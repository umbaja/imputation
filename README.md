# Genome Normalizer

Systém na spracovanie genotypových dát z rôznych čipových platforiem.
Jednu vzorku prevedie do spoločného normalizovaného genotypu (GRCh37, plus
vlákno, stĺpce `rsid/chromosome/position/genotype`) a vytvorí
dva oddelené výstupy:

1. **sibling/IBD dataset** — iba priamo namerané autosomálne genotypy;
2. **predispozičný panel** — 561 cieľových SNP, pričom chýbajúce lokusy možno
   cielene imputovať pomocou Beagle a referencie 1000 Genomes.

Genetické dáta možno spracovať lokálne alebo v heslom chránenej Railway službe.

Verejná inštancia: <https://genome-normalizer-production.up.railway.app/>

Bezpečné zdieľané demo používa iba syntetické dáta:

```text
používateľ: demo
heslo:      genome-demo-2026
```

DEMO účet nepovoľuje vlastné uploady ani Beagle imputáciu. Po kliknutí na
„Spustiť demo analýzu“ vykoná reálnu normalizáciu syntetického MyHeritage
súboru, panelovú kontrolu a ponúkne výsledky na stiahnutie.

Obmedzený TEST účet umožňuje spracovať vlastné dáta v ostrej pipeline:

```text
používateľ: test
heslo:      odovzdáva sa oprávneným testerom mimo verejného repozitára
limit:      5 genotypov denne, spoločne pre celý účet
```

Kvóta sa obnovuje o polnoci v časovom pásme Europe/Bratislava a jej stav
pretrváva na Railway volume aj po reštarte. Účet používajte iba na dáta, ktoré
máte oprávnenie spracúvať.

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

## Normalizovaný profil a kompatibilita v5

Bez parametra `--v5-template` program vytvorí normalizovaný genotyp (GRCh37,
plus vlákno, štyri stĺpce) a zachová všetky namerané markery. Táto otvorená
tabuľková schéma je kompatibilná s bežnými downstream nástrojmi. Repozitár
neobsahuje úplný licenčne overený manifest platformy v5.

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
popísané v [REFERENCE_DATA.md](REFERENCE_DATA.md). Úplná technická správa
GitHub/Railway prevádzky je v [TECHNICAL_REPORT.md](TECHNICAL_REPORT.md) a
v graficky spracovanej
[PDF verzii](output/pdf/Genome_Normalizer_Technicka_sprava_GitHub_Railway.pdf).
