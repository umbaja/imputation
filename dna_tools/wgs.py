"""
wgs.py — spracovanie WGS/WES dat (GATK VCF + BAM) do 23andMe formatu.

Preco samostatny modul: VCF zo sekvenovania NIE JE "dalsi cipovy format".
Lisi sa v troch podstatnych veciach:

  1. BUILD. Vacsina moderných pipeline bezi na hg38/GRCh38, zatial co cela tato
     appka (23andMe referencia, 1000G bref3, human_g1k_v37.fasta) je GRCh37/b37.
     Pozicie sa preto vobec netrafia -> tu robime liftover hg38 -> GRCh37.

  2. ZIADNE rsID. GATK HaplotypeCaller bez dbSNP anotacie necha ID stlpec '.',
     takze rsID musime doplnit podla POZICIE z 23andMe referencie.

  3. LEN VARIANTNE MIESTA. Homozygotne-referencne pozicie v subore nie su.
     Rozdiel medzi "je to ref/ref" a "nevieme, nebolo sekvenovane" sa da urobit
     LEN z pokrytia -> citame hlbku z BAM (samtools depth). Kde je hlbka >= MIN_DP
     a nie je tam variant, zapiseme homozygot referencie z FASTA; kde nie je
     pokrytie, ostava '--' (a moze ist do imputacie).

Výstupom je plnohodnotný normalizovaný genotyp nad šablónou pozícií z
data/id_position_reference_v3v4v5.tsv.gz, takze zvysok pipeline
(panel-check, cielena imputacia, davka) funguje bez zmeny.

Potrebuje: pyliftover (pip), chain subory z UCSC, samtools na PATH.
Vsetko doinstaluje/stiahne `bash setup_wgs.sh`.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Optional, Set, Tuple

ROOT = Path(__file__).resolve().parent.parent
CHAIN_DIR = ROOT / "ref" / "chains"
CHAIN_38TO19 = CHAIN_DIR / "hg38ToHg19.over.chain.gz"
CHAIN_19TO38 = CHAIN_DIR / "hg19ToHg38.over.chain.gz"
TEMPLATE_CACHE = ROOT / "data" / "template_hg38.tsv.gz"

_ACGT = set("ACGT")
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

# dlzka chr1 jednoznacne odlisi build
_CHR1_LEN = {248956422: "hg38", 249250621: "b37"}


# --- drobnosti ----------------------------------------------------------------

def _open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") \
        if str(path).endswith(".gz") else open(path, "r", encoding="utf-8", errors="replace")


def rc(seq: str) -> str:
    return "".join(_COMP.get(b, "N") for b in reversed(seq))


def to_ucsc(chrom: str) -> str:
    """'1' -> 'chr1', 'MT' -> 'chrM'."""
    c = str(chrom).strip()
    if c.startswith("chr"):
        return "chrM" if c == "chrMT" else c
    if c in ("MT", "M"):
        return "chrM"
    return "chr" + c


def from_ucsc(chrom: str) -> Optional[str]:
    """'chr1' -> '1', 'chrM' -> 'MT'. Ne-primarne kontigy (chr1_KI2707..) -> None."""
    c = str(chrom).strip()
    if c.startswith("chr"):
        c = c[3:]
    if "_" in c or c in ("EBV", "Un"):
        return None
    if c in ("M", "MT"):
        return "MT"
    if c in ("X", "Y") or (c.isdigit() and 1 <= int(c) <= 22):
        return c
    return None


# --- FASTA random access (cez .fai, bez pysam) --------------------------------

class FastaIndex:
    """Nacita <fasta>.fai a vie precitat jednu bazu na (chrom, 1-based pos).

    Zamerne bez pysam — .fai staci a je to par riadkov. Drzi otvoreny file
    handle a seekuje; pre ~1M dopytov je to radovo sekundy.
    """

    def __init__(self, fasta_path: str):
        self.path = str(fasta_path)
        fai = self.path + ".fai"
        if not Path(fai).exists():
            raise FileNotFoundError(
                f"Chyba index {fai}. Vytvor ho: samtools faidx {self.path}")
        self.idx: Dict[str, Tuple[int, int, int, int]] = {}
        with open(fai, "r", encoding="utf-8") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                # name, length, offset, linebases, linewidth
                self.idx[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
        self._fh = open(self.path, "rb")

    def _key(self, chrom: str) -> Optional[str]:
        c = str(chrom)
        for cand in (c, c.replace("chr", ""), "chr" + c.replace("chr", ""),
                     "MT" if c in ("chrM", "M") else c):
            if cand in self.idx:
                return cand
        return None

    def base(self, chrom: str, pos: int) -> Optional[str]:
        """1-based pozicia -> velke pismeno bazy, alebo None."""
        k = self._key(chrom)
        if k is None:
            return None
        length, offset, linebases, linewidth = self.idx[k]
        if pos < 1 or pos > length:
            return None
        p = pos - 1
        byte = offset + (p // linebases) * linewidth + (p % linebases)
        self._fh.seek(byte)
        b = self._fh.read(1).decode("ascii", "replace").upper()
        return b if b in _ACGT else (b if b else None)

    def close(self):
        try:
            self._fh.close()
        except Exception:  # noqa: BLE001
            pass


# --- detekcia buildu ----------------------------------------------------------

def detect_vcf_build(path: str) -> Tuple[str, str]:
    """Vrati ('hg38'|'b37'|'unknown', vysvetlenie) na zaklade hlavicky VCF."""
    ref_line = None
    with _open_text(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##reference="):
                ref_line = line.strip()
            if line.startswith("##contig=") and ("ID=chr1," in line or "ID=1," in line):
                for part in line.strip().rstrip(">").split(","):
                    if part.startswith("length="):
                        try:
                            ln = int(part.split("=")[1])
                        except ValueError:
                            continue
                        if ln in _CHR1_LEN:
                            return (_CHR1_LEN[ln],
                                    f"dĺžka chr1 = {ln} v hlavičke ##contig")
            if line.startswith("#CHROM"):
                break
    if ref_line:
        low = ref_line.lower()
        if "hg38" in low or "grch38" in low:
            return "hg38", f"hlavička {ref_line[:80]}"
        if "hg19" in low or "grch37" in low or "g1k_v37" in low or "b37" in low:
            return "b37", f"hlavička {ref_line[:80]}"
    return "unknown", "build sa z hlavičky VCF nedá určiť"


# --- liftover -----------------------------------------------------------------

def load_lifter(chain_path: Path):
    try:
        from pyliftover import LiftOver
    except ImportError as e:  # noqa: BLE001
        raise RuntimeError(
            "Chýba balík pyliftover. Nainštaluj ho: `bash setup_wgs.sh` "
            "(alebo pip install pyliftover v .venv-wsl)") from e
    if not Path(chain_path).exists():
        raise RuntimeError(
            f"Chýba chain súbor {chain_path}. Stiahni ho: `bash setup_wgs.sh`")
    return LiftOver(str(chain_path))


def lift_point(lifter, chrom_ucsc: str, pos1: int):
    """1-based pozicia -> (chrom_ucsc, pos1, strand) alebo None (nejednoznacne/nemapovatelne)."""
    res = lifter.convert_coordinate(chrom_ucsc, pos1 - 1)
    if not res:
        return None
    c, p0, strand = res[0][0], res[0][1], res[0][2]
    return c, p0 + 1, strand


# --- sablona 23andMe pozicii ---------------------------------------------------

def load_template(id_reference: str) -> list:
    """[(rsid, chrom37, pos37)] z 23andMe id_position_reference."""
    out = []
    with _open_text(id_reference) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h.strip().lower(): i for i, h in enumerate(header)}
        i_id = idx.get("id", idx.get("rsid", 0))
        i_chr = idx.get("chromosome", idx.get("chr", 1))
        i_pos = idx.get("position", idx.get("pos", 2))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_id, i_chr, i_pos):
                continue
            try:
                pos = int(float(f[i_pos]))
            except (ValueError, IndexError):
                continue
            out.append((f[i_id].strip(), f[i_chr].strip(), pos))
    return out


def build_template_hg38(id_reference: str, cache: Optional[Path] = None,
                        chain: Optional[Path] = None, progress=None) -> Path:
    """Raz prelozi vsetky sablonove pozicie GRCh37 -> hg38 a ulozi do cache."""
    cache = Path(cache or TEMPLATE_CACHE)      # rozlisenie az za behu (kvoli testom)
    chain = Path(chain or CHAIN_19TO38)
    if cache.exists():
        return cache
    lifter = load_lifter(chain)
    rows = load_template(id_reference)
    cache.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_suffix(".tmp.gz")
    n_ok = 0
    with gzip.open(tmp, "wt", encoding="utf-8", newline="\n") as fh:
        fh.write("rsid\tchrom37\tpos37\tchrom38\tpos38\n")
        for i, (rsid, c37, p37) in enumerate(rows):
            m = lift_point(lifter, to_ucsc(c37), p37)
            if m is None:
                continue
            c38, p38, _ = m
            if from_ucsc(c38) is None:          # spadlo na alt/random kontig
                continue
            fh.write(f"{rsid}\t{c37}\t{p37}\t{c38}\t{p38}\n")
            n_ok += 1
            if progress and i % 50000 == 0:
                progress(f"liftover šablóny: {i}/{len(rows)}")
    tmp.replace(cache)
    if progress:
        progress(f"šablóna hg38 hotová: {n_ok}/{len(rows)} pozícií -> {cache}")
    return cache


def load_template_hg38(cache: Optional[Path] = None) -> list:
    """[(rsid, chrom37, pos37, chrom38, pos38)]"""
    out = []
    with _open_text(Path(cache or TEMPLATE_CACHE)) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            out.append((f[0], f[1], int(f[2]), f[3], int(f[4])))
    return out


# --- VCF -> genotypy v GRCh37 --------------------------------------------------

@dataclass
class VcfStats:
    total: int = 0
    kept: int = 0
    skipped_nonsnv: int = 0
    skipped_filter: int = 0
    skipped_nogt: int = 0
    skipped_nolift: int = 0
    skipped_refmismatch: int = 0
    strand_flips: int = 0


def vcf_genotypes_b37(vcf_path: str, fai: FastaIndex, lifter=None,
                      pass_only: bool = True):
    """Precita VCF a vrati (genotypy, podozrive, statistika).

      genotypy  : {(chrom37, pos37): 'AG'} — spolahlive volania
      podozrive : {(chrom37, pos37)} — miesta, kde variant BOL, ale sme ho zahodili
                  (nepresiel filtrom / REF nesedi / chyba GT). Tieto pozicie sa
                  NESMU neskor vyhlasit za ref/ref len preto, ze tam nie je
                  genotyp — vieme, ze tam nieco je. Dostanu '--'.

    lifter=None znamena, ze VCF uz je v GRCh37 (neliftuje sa).
    Alely su vzdy na PLUS vlakne GRCh37 (pri strand flipe sa komplementuju).
    """
    gts: Dict[Tuple[str, int], str] = {}
    suspect: Set[Tuple[str, int]] = set()
    st = VcfStats()
    with _open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            st.total += 1
            chrom, pos_s, ref, alt, filt = f[0], f[1], f[3].upper(), f[4].upper(), f[6]
            if len(ref) != 1 or ref not in _ACGT:
                st.skipped_nonsnv += 1
                continue
            alts = alt.split(",")
            if any(len(a) != 1 or a not in _ACGT for a in alts):
                st.skipped_nonsnv += 1
                continue
            try:
                pos = int(pos_s)
            except ValueError:
                continue

            # --- liftover na GRCh37 (robime ho aj pre zahodene, aby sme vedeli
            #     oznacit podozrive pozicie) ---
            strand = "+"
            if lifter is None:
                c37, p37 = (from_ucsc(chrom) or str(chrom)), pos
            else:
                m = lift_point(lifter, to_ucsc(chrom), pos)
                if m is None:
                    st.skipped_nolift += 1
                    continue
                c38out, p37, strand = m
                c37 = from_ucsc(c38out)
                if c37 is None:
                    st.skipped_nolift += 1
                    continue

            # --- kontrola REF oproti GRCh37 FASTA ---
            ref37 = _COMP.get(ref, "N") if strand == "-" else ref
            base = fai.base(c37, p37)
            if base is None or base != ref37:
                # REF sa medzi buildmi lisi (alebo zla pozicia) -> zahodit.
                # ZAMERNE tu NIE JE "skus komplement" fallback: pri A/T a C/G
                # markeroch by ticho prijal nespravne alely. Chain uz strand riesi.
                st.skipped_refmismatch += 1
                suspect.add((c37, p37))
                continue

            if pass_only and filt not in ("PASS", ".", ""):
                st.skipped_filter += 1
                suspect.add((c37, p37))
                continue

            fmt = f[8].split(":")
            if "GT" not in fmt:
                st.skipped_nogt += 1
                suspect.add((c37, p37))
                continue
            gt = f[9].split(":")[fmt.index("GT")].replace("|", "/")
            alleles = [ref] + alts
            try:
                ix = [int(x) for x in gt.split("/") if x != "."]
                a = [alleles[i] for i in ix]
            except (ValueError, IndexError):
                a = []
            if not a:
                st.skipped_nogt += 1
                suspect.add((c37, p37))
                continue
            if len(a) == 1:
                a = [a[0], a[0]]          # haploid (chrX/Y u muza) -> hom

            if strand == "-":
                a = [_COMP.get(x, "N") for x in a]
                st.strand_flips += 1

            gts[(c37, p37)] = "".join(sorted(a))
            st.kept += 1
    return gts, suspect, st


# --- pokrytie z BAM ------------------------------------------------------------

def covered_positions(bam_path: str, positions38: Iterable[Tuple[str, int]],
                      min_dp: int = 8, min_mapq: int = 10, min_bq: int = 10,
                      progress=None) -> Set[Tuple[str, int]]:
    """`samtools depth` nad zoznamom hg38 pozicii -> mnozina tych s hlbkou >= min_dp."""
    if not shutil.which("samtools"):
        raise RuntimeError("samtools nie je na PATH — bez neho sa pokrytie z BAM nedá zistiť.")
    if not Path(bam_path).exists():
        raise RuntimeError(f"BAM súbor neexistuje: {bam_path}")
    pos_list = sorted(set(positions38), key=lambda t: (t[0], t[1]))
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False,
                                     encoding="utf-8", newline="\n") as bf:
        for c, p in pos_list:
            bf.write(f"{c}\t{p-1}\t{p}\n")
        bed = bf.name
    covered: Set[Tuple[str, int]] = set()
    try:
        cmd = ["samtools", "depth", "-a", "-b", bed,
               "-Q", str(min_mapq), "-q", str(min_bq), str(bam_path)]
        if progress:
            progress("čítam pokrytie z BAM: " + " ".join(cmd[:4]) + " …")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True, bufsize=1 << 20)
        n = 0
        for line in proc.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                if int(f[2]) >= min_dp:
                    covered.add((f[0], int(f[1])))
            except ValueError:
                continue
            n += 1
            if progress and n % 200000 == 0:
                progress(f"pokrytie: {n} pozícií prečítaných, {len(covered)} nad prahom")
        err = proc.stderr.read()
        proc.wait()
        if proc.returncode != 0:
            raise RuntimeError(f"samtools depth zlyhal ({proc.returncode}): {err.strip()[:300]}")
    finally:
        try:
            os.unlink(bed)
        except OSError:
            pass
    return covered


# --- orchestrator --------------------------------------------------------------

@dataclass
class WesResult:
    out_path: str = ""
    n_template: int = 0
    n_variant: int = 0
    n_homref: int = 0
    n_nocall: int = 0
    build: str = ""
    build_reason: str = ""
    bam_used: bool = False
    vcf_stats: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)


def wes_vcf_to_23andme(vcf_path: str, out_path: str, fasta_b37: str,
                       id_reference: str, bam_path: Optional[str] = None,
                       min_dp: int = 8, sample_name: str = "SAMPLE",
                       progress=None) -> WesResult:
    """Hlavná funkcia: (hg38) GATK VCF [+ BAM] -> normalizovaný genotyp."""
    def p(msg):
        if progress:
            progress(msg)

    res = WesResult()
    build, why = detect_vcf_build(vcf_path)
    res.build, res.build_reason = build, why
    p(f"build VCF: {build} ({why})")

    lifter = None
    if build == "hg38":
        p("načítavam chain hg38 → GRCh37…")
        lifter = load_lifter(CHAIN_38TO19)
    elif build == "unknown":
        res.warnings.append(
            "Build sa nepodarilo určiť z hlavičky — predpokladám GRCh37. "
            "Ak sú dáta v hg38, výsledok bude nezmyselný.")

    p("pripravujem šablónu 23andMe pozícií…")
    if build == "hg38":
        build_template_hg38(id_reference, progress=progress)
        tpl = load_template_hg38()
    else:
        tpl = [(rs, c, pos, to_ucsc(c), pos) for rs, c, pos in load_template(id_reference)]
    res.n_template = len(tpl)
    p(f"šablóna: {len(tpl)} pozícií")

    fai = FastaIndex(fasta_b37)
    try:
        p("čítam VCF…")
        gts, suspect, st = vcf_genotypes_b37(vcf_path, fai, lifter=lifter)
        res.vcf_stats = st.__dict__.copy()
        p(f"VCF: {st.kept} použiteľných variantov z {st.total} "
          f"(nemapované {st.skipped_nolift}, REF nesedí {st.skipped_refmismatch})")

        covered: Set[Tuple[str, int]] = set()
        if bam_path:
            covered = covered_positions(bam_path, [(c38, p38) for _, _, _, c38, p38 in tpl],
                                        min_dp=min_dp, progress=progress)
            res.bam_used = True
            p(f"BAM: {len(covered)} šablónových pozícií má hĺbku ≥ {min_dp}×")
        else:
            res.warnings.append(
                "BAM nebol zadaný — pozície bez variantu nevieme odlíšiť od nesekvenovaných, "
                "preto ostávajú ako no-call ('--') a pôjdu do imputácie.")

        p("zapisujem normalizovaný genotyp…")
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("# This data file generated by dna_tools (23andMe-compatible export).\n#\n")
            fh.write("# Zdroj: WGS/WES VCF%s.\n" % (" + BAM (pokrytie)" if bam_path else ""))
            fh.write("# Reference human assembly build 37 (GRCh37).\n")
            fh.write("# Build vstupneho VCF: %s (%s)\n" % (build, why))
            fh.write("# Sample: %s\n#\n" % sample_name)
            fh.write("# rsid\tchromosome\tposition\tgenotype\n")
            for rsid, c37, p37, c38, p38 in tpl:
                gt = gts.get((c37, p37))
                if gt:
                    res.n_variant += 1
                elif (c37, p37) in suspect:
                    # variant tam JE, len sme mu neverili -> radsej no-call ako ref/ref
                    gt = "--"
                    res.n_nocall += 1
                elif covered and (c38, p38) in covered:
                    b = fai.base(c37, p37)
                    if b in _ACGT:
                        gt = b + b
                        res.n_homref += 1
                    else:
                        gt = "--"
                        res.n_nocall += 1
                else:
                    gt = "--"
                    res.n_nocall += 1
                fh.write(f"{rsid}\t{c37}\t{p37}\t{gt}\n")
    finally:
        fai.close()

    res.out_path = str(out_path)
    p(f"hotovo: {res.n_variant} variantov, {res.n_homref} ref/ref, {res.n_nocall} no-call")
    return res


def wgs_status() -> dict:
    """Co je / nie je pripravene pre WGS vetvu (pre /api/status)."""
    try:
        import pyliftover  # noqa: F401
        have_lift = True
    except ImportError:
        have_lift = False
    return {
        "pyliftover": have_lift,
        "samtools": bool(shutil.which("samtools")),
        "chain_38to19": CHAIN_38TO19.exists(),
        "chain_19to38": CHAIN_19TO38.exists(),
        "template_hg38": TEMPLATE_CACHE.exists(),
        "wgs_ready": have_lift and CHAIN_38TO19.exists() and CHAIN_19TO38.exists(),
    }
