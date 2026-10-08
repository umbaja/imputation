# Genome Normalizer

Lokálny systém na spracovanie genotypových dát z rôznych čipových platforiem.
Jednu vzorku normalizuje do spoločnej reprezentácie GRCh37/23andMe a vytvorí
dva oddelené výstupy:

1. **sibling/IBD dataset** — iba priamo namerané autosomálne genotypy;
2. **predispozičný panel** — 561 cieľových SNP, pričom chýbajúce lokusy možno
   cielene imputovať pomocou Beagle a referencie 1000 Genomes.

Genetické dáta sa štandardne spracúvajú iba na lokálnom počítači.

## Podporované vstupy

- 23andMe raw data
- AncestryDNA
- MyHeritage a FamilyTreeDNA
- Illumina GenomeStudio Final Report
- VCF
- PLINK PED/MAP
- WGS/WES VCF s voliteľným BAM vstupom

Illumina dáta musia pre spoločný analytický výstup obsahovať `Allele1 - Plus`
a `Allele2 - Plus`. TOP/BOT alebo Forward export bez manifestovej konverzie je
zámerne odmietnutý, pretože by mohol obrátiť alely.

## Najrýchlejšie spustenie

Na Windowse sa odporúča WSL2:

```bash
git clone https://github.com/umbaja/imputation.git
cd imputation
bash install_tools_root.sh       # raz, so sudo/root oprávnením
bash setup_reference.sh          # raz; stiahne a SHA-256 overí referencie
bash start_app.sh
```

Webové rozhranie bude na <http://127.0.0.1:8000>.

Postup pre heslom chránené verejné nasadenie na Railway je v
[RAILWAY_DEPLOY.md](RAILWAY_DEPLOY.md). Referenčné dáta ostávajú mimo GitHubu
a pri prvom štarte sa stiahnu na pripojený persistentný volume.

Bez imputácie stačí Python 3.9+:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
uvicorn app.main:app --host 127.0.0.1 --port 8000
```

## Jeden vstup, dva analytické výstupy

```bash
genome-normalizer prepare-analysis VZORKA.txt --output-dir out/vzorka
```

So zapnutou cielenou imputáciou predispozičného panela:

```bash
export BEAGLE_JAR="$PWD/ref/beagle.jar"
export REF_DIR="$PWD/ref/b37.bref3"
export MAP_DIR="$PWD/ref/maps"
export FASTA="$PWD/ref/human_g1k_v37.fasta"

genome-normalizer prepare-analysis VZORKA.txt \
  --output-dir out/vzorka \
  --impute \
  --min-gp 0.9
```

Každý beh vytvorí normalizovaný súbor, sibling dataset, stav 561-SNP panela,
výsledný predispozičný panel a JSON QC report. Imputované hodnoty sa nikdy
nepridávajú do sibling/IBD výstupu.

## Profil 23andMe v5

Bez parametra `--v5-template` program vytvorí 23andMe-v5-kompatibilnú schému
(GRCh37, štyri stĺpce) a zachová všetky namerané markery. Repozitár neobsahuje
úplný licenčne overený manifest 23andMe v5.

Ak je dostupný presný v5 markerový súbor, možno vytvoriť pevnú kostru:

```bash
genome-normalizer prepare-analysis VZORKA.txt \
  --output-dir out/vzorka \
  --v5-template /cesta/23andme-v5-template.txt
```

Nedostupné v5 markery ostanú ako no-call. Predispozičný panel sa stále
vyhodnocuje proti plnému normalizovanému vstupu, aby sa nestratili priamo
namerané panelové SNP mimo v5.

## Docker

Referencie sa nezapisujú do obrazu; pripájajú sa ako volume:

```bash
docker build -t genome-normalizer .
docker run --rm -p 8000:8000 \
  -e PUBLIC_MODE=0 \
  -v "$PWD/ref:/app/ref:ro" \
  -v "$PWD/out:/app/out" \
  genome-normalizer
```

## Dôležité obmedzenia

- Súradnicový základ je GRCh37/hg19. Hg38 vstupy musia prejsť validovaným
  liftoverom.
- Pre sibling analýzu sa používajú iba reálne namerané genotypy. Ak majú dve
  platformy menej než 10 000 spoločných autosomálnych SNP, cielená imputácia
  561 panelových lokusov tento problém nevyrieši.
- Predispozičné výsledky sú výskumné/informatívne. Klinicky významný nález sa
  musí potvrdiť nezávislým laboratórnym vyšetrením.
- Referencie Beagle/1000G majú približne 13 GB a spravuje ich
  `setup_reference.sh`; do Git repozitára nepatria. Presné verzie, zdroje,
  veľkosti a SHA-256 sú pripnuté v `reference_manifest.tsv`. Už stiahnutú
  referenciu možno bez siete skontrolovať cez
  `bash setup_reference.sh --verify-only`.

Podrobný návrh komponentov a dátového toku je v
[ARCHITECTURE.md](ARCHITECTURE.md). Bezpečnostné pravidlá sú v
[SECURITY.md](SECURITY.md). Pôvod a správa veľkých externých súborov sú
popísané v [REFERENCE_DATA.md](REFERENCE_DATA.md).
