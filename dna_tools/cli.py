"""
cli.py
======

Príkazový riadok pre dna_tools.

Príklady:
    python -m dna_tools detect  vstup.txt
    python -m dna_tools convert vstup.txt -o vystup_23andme.txt
    python -m dna_tools convert illumina_report.txt -o out.txt --drop-no-call
    python -m dna_tools compare  suborA.txt suborB.csv --key pos --mismatch-csv nezhody.csv
"""

from __future__ import annotations

import argparse
import sys

from .detect import detect_format_verbose
from .parsers import parse_file
from .convert import write_23andme
from .compare import compare_files
from .impute import build_scaffold
from .panel import check_panel
from .vcfio import write_vcf_from_sample
from .targeted_impute import build_windows, write_regions_bed, extract_targets_from_vcf, merge_panel
from .convert import write_23andme as _write_23
from .analysis_prep import prepare_for_analyses


def _cmd_detect(args) -> int:
    fmt, reason = detect_format_verbose(args.input)
    print(f"Súbor : {args.input}")
    print(f"Formát: {fmt.value}")
    print(f"Dôvod : {reason}")
    return 0 if fmt.value != "unknown" else 2


def _cmd_convert(args) -> int:
    data = parse_file(args.input, map_path=args.map)
    print(data.summary(), file=sys.stderr)
    for w in data.warnings:
        print(f"  ⚠ {w}", file=sys.stderr)
    out = write_23andme(
        data, args.output,
        drop_no_call=args.drop_no_call,
        keep_non_rs=not args.only_rs,
        sort_alleles=not args.no_sort_alleles,
    )
    print(f"Zapísané do: {out}")
    return 0


def _cmd_compare(args) -> int:
    result = compare_files(
        args.file_a, args.file_b,
        key=args.key,
        mismatch_csv=args.mismatch_csv,
    )
    print(result.summary())
    if args.mismatch_csv:
        print(f"\nNezhody zapísané do: {args.mismatch_csv}")
    return 0


def _cmd_scaffold(args) -> int:
    gd, stats = build_scaffold(args.sample, args.template,
                               sort_alleles=not args.no_sort_alleles)
    for w in gd.warnings:
        print(f"  ℹ {w}", file=sys.stderr)
    _write_23(gd, args.output, sort_alleles=False)  # uz zoradene v scaffolde
    print(f"Zapisane do: {args.output}")
    print("Pokrytie: %d/%d vyplnenych (%.1f%%), %d chyba (no-call)."
          % (stats["filled_from_sample"], stats["template_markers"],
             100.0*stats["fill_rate"], stats["missing_no_call"]))
    return 0


def _cmd_panel_check(args) -> int:
    res = check_panel(args.sample, args.panel)
    print(res.summary())
    if args.out_missing:
        miss = res.table[res.table["status"] != "called"]
        miss.to_csv(args.out_missing, index=False)
        print(f"\nChybajuce/no-call lokusy zapisane do: {args.out_missing}")
    if args.out_all:
        res.table.to_csv(args.out_all, index=False)
        print(f"Kompletny prehlad zapisany do: {args.out_all}")
    return 0


def _cmd_impute_targets(args) -> int:
    import csv as _csv
    targets=[(r["chromosome"], int(float(r["position"])))
             for r in _csv.DictReader(open(args.targets)) if r.get("position")]
    regs=build_windows(targets, flank=args.flank)
    total=sum(r.end-r.start for r in regs)
    print(f"Cielov: {len(targets)} | regionov (+-{args.flank} bp): {len(regs)} | pokryte: {total/1e6:.1f} Mb")
    if args.regions_bed:
        write_regions_bed(regs, args.regions_bed)
        print(f"Regiony (BED): {args.regions_bed}")
    if args.imputed_vcf:
        tmap={(c, p): rs for (c, p), rs in
              ((( r["chromosome"], int(float(r["position"]))), r["rsid"])
               for r in _csv.DictReader(open(args.targets)) if r.get("position"))}
        imp=extract_targets_from_vcf(args.imputed_vcf, tmap, min_gp=args.min_gp)
        print(f"Extrahovanych cielov z VCF: {len(imp)} | neiste (max_gp<{args.min_gp}): {int(imp['low_conf'].sum())}")
        if args.panel_status and args.output:
            final=merge_panel(args.panel_status, imputed=imp)
            final["position"]=final["position"].astype(float).astype("Int64")
            final.to_csv(args.output, index=False)
            print("Zdroje:", final["source"].value_counts().to_dict())
            print(f"Finalny panel: {args.output}")
    else:
        print("Dalej: spusti Beagle cez targeted_impute.sh (potrebuje ref. panel),")
        print("potom znova s --imputed-vcf ... --panel-status ... -o final.csv")
    return 0


def _cmd_to_vcf(args) -> int:
    from .vcfio import extract_sites_vcf
    chroms = args.chrom.split(",") if args.chrom else None
    sites = extract_sites_vcf(args.sites_vcf) if args.sites_vcf else None
    stats = write_vcf_from_sample(args.input, args.output, fasta_path=args.fasta,
                                  sites=sites, chroms=chroms)
    print("VCF: %s | zapisanych: %d | preskoc.: homref %d, ne-SNP %d, mismatch %d, bez-site %d"
          % (args.output, stats["written"], stats["skipped_homref"],
             stats["skipped_nonsnp"], stats["skipped_mismatch"], stats["skipped_nosite"]))
    return 0


def _cmd_prepare_analysis(args) -> int:
    result = prepare_for_analyses(
        sample_path=args.input,
        output_dir=args.output_dir,
        panel_path=args.panel,
        id_reference_path=args.id_reference,
        v5_template_path=args.v5_template,
        online_annotation=args.online_annotation,
        run_imputation=args.impute,
        min_gp=args.min_gp,
        flank=args.flank,
        beagle_jar=args.beagle_jar,
        ref_dir=args.ref_dir,
        map_dir=args.map_dir,
        fasta=args.fasta,
    )
    print(f"Normalizovaný dataset : {result.normalized_dataset}")
    print(f"Sibling dataset       : {result.sibling_dataset}")
    print(f"Predispozičný panel   : {result.predisposition_panel}")
    if result.predisposition_dataset:
        print(f"Dataset s imputáciou  : {result.predisposition_dataset}")
    print(f"QC report             : {result.report}")
    print(
        "Panel: %d meraných, %d imputovaných, %d odmietnutých pre nízku istotu, %d chýbajúcich"
        % (
            result.panel_measured,
            result.panel_imputed,
            result.panel_rejected_low_conf,
            result.panel_missing_after,
        )
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="dna_tools",
        description="Detekcia formátu, porovnanie pozícií a prevod genetických dát do 23andMe formátu.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    pd_ = sub.add_parser("detect", help="Rozpozná formát vstupného súboru.")
    pd_.add_argument("input")
    pd_.set_defaults(func=_cmd_detect)

    pc = sub.add_parser("convert", help="Prevedie vstup do 23andMe formátu.")
    pc.add_argument("input")
    pc.add_argument("-o", "--output", required=True, help="Cesta k výstupnému súboru.")
    pc.add_argument("--map", help="PLINK .map súbor (ak vstup je .ped).")
    pc.add_argument("--drop-no-call", action="store_true", help="Vynechať no-call ('--') pozície.")
    pc.add_argument("--only-rs", action="store_true", help="Ponechať len varianty s rsID.")
    pc.add_argument("--no-sort-alleles", action="store_true",
                    help="Nezoraďovať alely abecedne (predvolene sa zoraďujú kvôli 23andMe konvencii).")
    pc.set_defaults(func=_cmd_convert)

    pcmp = sub.add_parser("compare", help="Porovná zhodné pozície/genotypy dvoch súborov.")
    pcmp.add_argument("file_a")
    pcmp.add_argument("file_b")
    pcmp.add_argument("--key", choices=["rsid", "pos"], default="rsid",
                      help="Kľúč porovnania: rsID alebo chr:pos.")
    pcmp.add_argument("--mismatch-csv", help="Kam zapísať nezhody (CSV).")
    pcmp.set_defaults(func=_cmd_compare)

    psc = sub.add_parser("scaffold",
        help="Namapuj vzorku na marker set sablony (napr. 23andMe v5); chybajuce -> no-call.")
    psc.add_argument("sample", help="Vstupna vzorka (napr. MyHeritage).")
    psc.add_argument("--template", required=True, help="Sablona marker setu (napr. 23andMe v5 subor).")
    psc.add_argument("-o", "--output", required=True, help="Vystupny 23andMe subor.")
    psc.add_argument("--no-sort-alleles", action="store_true", help="Nezoradovat alely.")
    psc.set_defaults(func=_cmd_scaffold)

    ppc = sub.add_parser("panel-check",
        help="Overi pokrytie cieloveho panelu lokusov (~1000 pozicii) v subore.")
    ppc.add_argument("sample", help="Genotypovy subor (napr. prevedeny 23andMe).")
    ppc.add_argument("--panel", required=True, help="Panel: zoznam rsID alebo CSV/TSV (rsid/chrom/pos/effect_allele).")
    ppc.add_argument("--out-missing", help="CSV s chybajucimi/no-call lokusmi (kandidati na imputaciu).")
    ppc.add_argument("--out-all", help="CSV s prehladom vsetkych lokusov a ich stavom.")
    ppc.set_defaults(func=_cmd_panel_check)

    pit = sub.add_parser("impute-targets",
        help="Cielena imputacia: vypocet okien okolo cielov + zlucenie imputovanych genotypov.")
    pit.add_argument("--targets", required=True, help="CSV s cielmi (rsid,chromosome,position).")
    pit.add_argument("--flank", type=int, default=250000, help="Polsirka okna v bp (default 250000).")
    pit.add_argument("--regions-bed", help="Kam zapisat regiony (BED).")
    pit.add_argument("--imputed-vcf", help="Imputovany VCF z Beagle (na extrakciu cielov).")
    pit.add_argument("--panel-status", help="panel-check --out-all CSV (na zlucenie).")
    pit.add_argument("--min-gp", type=float, default=0.9, help="Prah max(GP) pre spolahlivost (default 0.9).")
    pit.add_argument("--dr2-min", type=float, default=0.3, help="(zastarale) Prah DR2 — pouzi --min-gp.")
    pit.add_argument("-o", "--output", help="Finalny zluceny panel CSV.")
    pit.set_defaults(func=_cmd_impute_targets)

    ptv = sub.add_parser("to-vcf",
        help="Prevedie 23andMe -> VCF. Rezim --fasta (informativne) alebo --sites-vcf (aj hom-ref).")
    ptv.add_argument("input", help="23andMe subor.")
    ptv.add_argument("--fasta", help="Referencna FASTA (rezim A: REF z FASTA).")
    ptv.add_argument("--sites-vcf", help="VCF s panelovymi sites (rezim B: aj hom-ref, viac kotiev).")
    ptv.add_argument("--chrom", help="Obmedz na chromozom(y), oddelene ciarkou.")
    ptv.add_argument("-o", "--output", required=True, help="Vystupny VCF.")
    ptv.set_defaults(func=_cmd_to_vcf)

    ppa = sub.add_parser(
        "prepare-analysis",
        help="Pripraví jednu vzorku pre sibling analýzu aj predispozičný panel.",
    )
    ppa.add_argument("input", help="Vstupná genotypová vzorka z podporovanej platformy.")
    ppa.add_argument("--output-dir", required=True, help="Priečinok pre všetky výstupy a QC report.")
    ppa.add_argument(
        "--panel",
        help="Predispozičný panel; default data/panel_full_561.csv.",
    )
    ppa.add_argument(
        "--id-reference",
        help="Lokálna GRCh37 referencia rsID/pozícií; default data/id_position_reference_v3v4v5.tsv.gz.",
    )
    ppa.add_argument(
        "--v5-template",
        help="Voliteľný presný 23andMe v5 markerový súbor pre sibling scaffold.",
    )
    ppa.add_argument(
        "--online-annotation",
        action="store_true",
        help="Dohľadať neznáme rsID cez Ensembl GRCh37.",
    )
    ppa.add_argument(
        "--impute",
        action="store_true",
        help="Spustiť Beagle iba pre chýbajúce SNP predispozičného panela.",
    )
    ppa.add_argument("--min-gp", type=float, default=0.9, help="Minimálne max(GP) pre zlúčenie imputácie.")
    ppa.add_argument("--flank", type=int, default=250000, help="Polovica Beagle okna v bp.")
    ppa.add_argument("--beagle-jar", help="Cesta k beagle.jar.")
    ppa.add_argument("--ref-dir", help="Priečinok 1000G b37 bref3.")
    ppa.add_argument("--map-dir", help="Priečinok genetických máp GRCh37.")
    ppa.add_argument("--fasta", help="Referenčná FASTA GRCh37.")
    ppa.set_defaults(func=_cmd_prepare_analysis)

    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
