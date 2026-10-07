# Samostatné fungovanie Genome Normalizer

## Cieľ

Repozitár funguje ako samostatná lokálna služba medzi surovými čipovými dátami
a dvoma konzumentmi:

```text
rôzne čipové platformy
        │
        ▼
detekcia formátu → GRCh37/plus normalizácia → QC a pôvod genotypu
        │
        ├──► sibling dataset (iba merané SNP) ──► PLINK/IBIS
        │
        └──► 561-SNP panel
                 ├── merané SNP
                 └── chýbajúce SNP ──► Beagle/1000G ──► predispozície
```

Jeden normalizovaný zdroj je spoločný, ale analytické výstupy sú oddelené.
Toto je dôležité: cielená panelová imputácia nesmie byť vydávaná za merané
dáta pri IBD analýze.

## Komponenty

- `dna_tools/` — detekcia formátu, parsovanie, normalizácia, panel-check,
  VCF konverzia, extrakcia Beagle výstupu a CLI;
- `dna_tools/analysis_prep.py` — spoločná orchestrace dvoch analytických vetiev;
- `targeted_impute.sh` — izolovaný Beagle beh po cieľových regiónoch;
- `app/` a `static/` — lokálne FastAPI webové rozhranie;
- `data/panel_full_561.csv` — verzovaný predispozičný panel;
- `data/id_position_reference_v3v4v5.tsv.gz` — malá lokálna mapa rsID/GRCh37;
- `ref/` — neverzované Beagle, 1000G, mapy a FASTA dáta.

## Prevádzkové režimy

### 1. Lokálne webové UI

Odporúčaný režim pre jednotlivé vzorky a používateľov bez príkazového riadka.
Služba počúva iba na `127.0.0.1`; surové genotypy sa neodosielajú na internet.

### 2. CLI a dávkové spracovanie

`genome-normalizer prepare-analysis` je stabilný integračný bod. Každá vzorka
má vlastný výstupný a pracovný priečinok, takže sa behy navzájom neprepisujú.
JSON report umožňuje nadväzujúcemu Sibling Analyzeru automaticky nájsť správny
výstup a overiť QC.

### 3. Docker

Kód je v obraze, veľké referencie a výsledky sú volumes. Kontajner je vhodný na
lokálny server alebo internú infraštruktúru. Verejný multi-user hosting sa bez
autentifikácie, šifrovania a pravidiel retencie genetických dát neodporúča.

## Integrácia so Sibling Analyzerom

Sibling Analyzer má dostať iba cestu `sibling_dataset` z JSON reportu. Pred
spustením PLINK/IBIS má overiť:

1. obe vzorky sú GRCh37;
2. obe prešli strand kontrolou;
3. majú aspoň 10 000 spoločných zavolaných autosomálnych SNP;
4. spoločné markery sa párujú podľa `chromosome:position`, nie názvu sondy.

Ak prienik nestačí, systém má výsledok označiť ako „nedostatočné dáta“. Nemá
automaticky používať 561 panelovo imputovaných SNP na IBD. Budúca voliteľná
vetva pre hustú kinship imputáciu musí mať vlastný validačný protokol, spoločný
markerový set, GP filter a porovnanie na známych príbuzenských pároch.

## Výstupy jedného behu

- `*.normalized_grch37_23andme.txt` — kanonická reprezentácia vstupu;
- `*.sibling_23andme_v5.txt` — iba pri dodanej presnej v5 šablóne;
- `*.predisposition_status.csv` — stav všetkých panelových SNP;
- `*.predisposition_missing.csv` — kandidáti na imputáciu;
- `*.predisposition_panel.csv` — merané a kvalitné imputované výsledky;
- `*.with_predisposition_imputation.txt` — rozšírený normalizovaný súbor;
- `*.preparation_report.json` — pôvod, počty, cesty a varovania.

## Odporúčaný ďalší vývoj

1. Dodať licenčne overený presný manifest 23andMe v5 alebo vlastný otvorený
   spoločný markerový profil.
2. Pridať manifestovú konverziu Illumina TOP/BOT na plus vlákno.
3. Zaviesť explicitnú detekciu buildu pre všetky array formáty a kontrolovaný
   hg38→GRCh37 liftover.
4. Napísať validačnú sadu s rovnakými vzorkami z viacerých platforiem.
5. Napojiť JSON report priamo na Sibling Analyzer namiesto manuálneho uploadu.
