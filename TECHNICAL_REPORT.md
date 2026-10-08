# Technická správa: Genome Normalizer

Stav dokumentu: 8. októbra 2026  
Verzia aplikácie: 2.3  
Repozitár: <https://github.com/umbaja/imputation>  
Produkcia: <https://genome-normalizer-production.up.railway.app/>

## 1. Účel systému

Genome Normalizer prijíma genotypové súbory z rôznych platforiem a vytvorí
jednotný normalizovaný genotyp v zostave GRCh37 na plus vlákne. Kanonický
textový výstup má štyri tabulátorom oddelené stĺpce:
`rsid`, `chromosome`, `position`, `genotype`. No-call hodnota je `--`.

Z jedného normalizovaného zdroja vznikajú dve striktne oddelené analytické
vetvy:

1. sibling/IBD dataset obsahuje iba priamo namerané genotypy;
2. predispozičný dataset používa vstavaný 561-SNP panel a môže doplniť chýbajúce
   panelové lokusy cielenou imputáciou Beagle/1000 Genomes.

Imputované lokusy sa nikdy nepridávajú do sibling/IBD datasetu. Toto pravidlo
bráni tomu, aby sa panelové odhady nesprávne použili ako meraný dôkaz
príbuznosti.

## 2. Podporované vstupy

- 23andMe raw data;
- AncestryDNA;
- MyHeritage a FamilyTreeDNA;
- Illumina GenomeStudio Final Report;
- VCF;
- PLINK PED/MAP;
- WGS/WES VCF s voliteľným BAM pri lokálnej prevádzke.

Formát sa deteguje z autoritatívnej hlavičky a štruktúry stĺpcov, nie iba z
názvu súboru. Illumina musí poskytovať alely na plus vlákne alebo mať
validovanú manifestovú konverziu. Neistý TOP/BOT export sa odmietne.

## 3. Dátový tok

```text
surový genotyp
  -> detekcia formátu
  -> parser platformy
  -> normalizácia chromozómu, pozície, alel a no-call hodnôt
  -> normalizovaný GRCh37/plus genotyp
       |-> sibling/IBD: iba merané lokusy
       `-> 561-SNP predispozičný panel
             |-> measured
             |-> no_call/absent
             `-> Beagle + 1000G -> imputed, iba ak max(GP) >= 0,90
```

Normalizácia zoraďuje chromozómy a pozície, stabilizuje diploidné alely a
zachováva pôvodné namerané volania. Panelová kontrola páruje prednostne podľa
rsID a následne podľa `chromosome:position`. Lokálna ID referencia môže doplniť
chýbajúcu súradnicu bez odoslania genotypu na internet.

Pri imputácii `targeted_impute.sh` vytvorí VCF, zoradí a indexuje ho, zoskupí
ciele do okien, pre príslušné chromozómy spustí Beagle, zlúči regióny a
extrahuje iba cieľové SNP. Predvolený prah je max(GP) 0,90. Namerané volanie má
vždy prednosť pred imputovaným.

## 4. Výstupy

- `*_normalized_genotype.txt` alebo `*.normalized_grch37.txt` – normalizovaný
  genotyp;
- `*_sibling_measured_normalized.txt` – sibling/IBD vstup bez imputácie;
- `*_predisposition_normalized.txt` – normalizovaný predispozičný výstup;
- `*_predisposition_status.csv` – stav každého SNP pred imputáciou;
- `*_predisposition_missing.csv` – lokusy navrhnuté na imputáciu;
- `*_predisposition_panel.csv` – merané, prijaté imputované a chýbajúce SNP;
- `*.preparation_report.json` – QC, počty, pôvod a upozornenia.

## 5. Softvérová architektúra

| Komponent | Zodpovednosť |
|---|---|
| `app/main.py` | FastAPI, middleware, jednotlivé API kroky, DEMO, downloady |
| `app/batch.py` | kompletná viac-kroková pipeline a progres dávky |
| `app/runtime_security.py` | autentifikácia, roly, limity, retencia, globálny job lock |
| `dna_tools/detect.py` | automatická detekcia vstupného formátu |
| `dna_tools/parsers.py` | parsery platforiem do `GenotypeData` |
| `dna_tools/convert.py` | normalizovaný štvorstĺpcový export a merge imputácie |
| `dna_tools/panel.py` | kontrola pokrytia 561-SNP alebo vlastného panela |
| `dna_tools/analysis_prep.py` | stabilný CLI integračný bod pre obe vetvy |
| `targeted_impute.sh` | cielená Beagle/1000G imputácia |
| `setup_reference.sh` | stiahnutie, veľkostná a SHA-256 kontrola referencií |
| `deploy/start.sh` | Railway bootstrap a štart Uvicornu na `$PORT` |
| `static/index.html` | samostatné webové UI bez frontend frameworku |

Pythonová časť používa FastAPI, Uvicorn, pandas, python-multipart a pyliftover.
Kontajner dopĺňa Java runtime, bcftools, bgzip/tabix, samtools, curl, unzip a
gzip.

## 6. API

| Metóda a cesta | Funkcia |
|---|---|
| `GET /healthz` | verejný Railway health check |
| `GET /api/status` | verzia, rola, referencia a pripravenosť nástrojov |
| `POST /api/detect` | detekcia vstupného formátu |
| `POST /api/convert` | vytvorenie normalizovaného genotypu |
| `POST /api/panel-check` | pokrytie panelových SNP |
| `POST /api/impute-plan` | výpočet cieľových okien |
| `POST /api/impute-run` | asynchrónne spustenie imputácie |
| `GET /api/impute-progress/{id}` | progres imputácie |
| `POST /api/merge-final` | vytvorenie predispozičného výstupu |
| `POST /api/batch-run` | kompletná pipeline pre 1–5 vzoriek |
| `GET /api/batch-progress/{id}` | progres dávky |
| `POST /api/batch-cancel/{id}` | požiadavka na zrušenie dávky |
| `POST /api/demo-run` | syntetická normalizácia a panelová kontrola |
| `GET /api/demo-download/{kind}` | pevne povolené DEMO súbory |
| `GET /api/download/{name}` | výstupy plného účtu z pracovného adresára |

Vo verejnom režime sú Swagger/OpenAPI a `/api/logs` vypnuté.

## 7. DEMO účet

Produkčný DEMO účet je určený na verejnú prezentáciu:

```text
používateľ: demo
heslo:      genome-demo-2026
```

Súbor `demo/demo_myheritage_raw.csv` je syntetický a nepatrí žiadnej osobe.
Kliknutie na „Spustiť demo analýzu“ používa reálny parser MyHeritage,
normalizačný export a kontrolu vstavaného 561-SNP panela. Používateľ vidí
náhľad, počty a môže stiahnuť vstup, normalizovaný genotyp, sibling dataset,
panel CSV a sanitizovanú JSON správu.

DEMO middleware povoľuje iba `/`, `/api/status`, `/api/demo-run`, DEMO
downloady a statické súbory. Vlastné uploady, všeobecné downloady, logy,
dávková pipeline a Beagle sú blokované odpoveďou HTTP 403. Beagle sa v deme
zámerne nespúšťa, aby zdieľaný účet nemohol spotrebovať výpočtové zdroje.

## 8. GitHub

Zdrojový repozitár je verejný `umbaja/imputation`; produkčná vetva je `main`.
Push do `main` spustí dve nezávislé udalosti:

1. GitHub Actions workflow `.github/workflows/tests.yml` nainštaluje balík na
   Python 3.11, spustí všetky `unittest` testy a overí syntax troch shell
   skriptov;
2. Railway GitHub integrácia zostaví nový Docker image a nasadí ho do
   produkčného prostredia.

Do Gitu sa neukladajú surové genetické dáta, pracovné výstupy ani veľké
referencie. `.gitignore` blokuje okrem iného `work/`, `out/`, `results/`,
`samples/`, `uploads/`, `ref/`, VCF, BAM, PLINK a odvodené indexy. Verzované sú
iba kód, malý 561-SNP panel, lokálna komprimovaná ID mapa, syntetické demo a
manifest referencií.

## 9. Railway produkcia

| Položka | Hodnota |
|---|---|
| projekt | `genome-normalizer` |
| environment | `production` |
| service | `genome-normalizer` |
| doména | `genome-normalizer-production.up.railway.app` |
| build | koreňový `Dockerfile`, Python 3.11 slim |
| start | `bash deploy/start.sh` |
| port | Railway premenná `$PORT`, aktuálne smerovaná na 8080 |
| volume | `genome-normalizer-volume` |
| mount | `/app/ref` |
| kapacita | 50 000 MB; využitie po bootstrapovaní približne 12 128 MB |
| repliky | 1 |
| región | Railway predvolený US East |

### Premenné prostredia

| Premenná | Produkčná hodnota alebo význam |
|---|---|
| `PUBLIC_MODE` | `1` |
| `APP_USERNAME` | vlastník; hodnota nie je dokumentovaná ako verejný údaj |
| `APP_PASSWORD` | tajný Railway secret |
| `DEMO_USERNAME` | `demo` |
| `DEMO_PASSWORD` | `genome-demo-2026` – zámerne verejné demo heslo |
| `BOOTSTRAP_REFERENCE` | `1` |
| `REFERENCE_DIR` | `/app/ref` |
| `BEAGLE_JAR` | `/app/ref/beagle.jar` |
| `REF_DIR` | `/app/ref/b37.bref3` |
| `MAP_DIR` | `/app/ref/maps` |
| `FASTA` | `/app/ref/human_g1k_v37.fasta` |
| `JOBS` | `1` |
| `XMX_JOB` | `3g` |
| `RUNTIME_RETENTION_HOURS` | `6` |
| `MAX_UPLOAD_MB` | `100` |
| `MAX_REQUEST_MB` | `250` |
| `MAX_BATCH_SAMPLES` | `5` |
| `RAILWAY_HEALTHCHECK_TIMEOUT_SEC` | `3600` |

`APP_PASSWORD` sa nesmie commitnúť, zapisovať do technickej správy ani
zverejniť v logu. Zmena vlastníckeho hesla sa robí iba cez Railway Variables a
vyvolá nové nasadenie.

## 10. Referenčné dáta a prvý štart

`reference_manifest.tsv` pripína URL, presnú veľkosť a SHA-256 každého súboru:

- Beagle JAR;
- genetické mapy GRCh37;
- 23 bref3 panelov pre chromozómy 1–22 a X z 1000 Genomes phase 3;
- komprimovanú a rozbalenú `human_g1k_v37.fasta`.

Pri prvom štarte `deploy/start.sh` overí, či je referencia kompletná. Ak nie,
spustí `setup_reference.sh`. Download používa resume, každý súbor sa pred
presunom z `.part` overí veľkosťou a SHA-256, mapy sa rozbalia a samtools
vytvorí FASTA index `.fai`. Aplikácia sa spustí až po úspešnej kontrole.

Volume je persistentný, preto ďalší deployment iba skontroluje prítomnosť
Beagle, FASTA, `.fai`, 23 bref3 panelov a 23 máp. Kód a `/app/work` sú na
ephemerálnom filesystéme kontajnera; referencie zostávajú na volume.

## 11. Bezpečnosť a životný cyklus dát

- verejný režim zlyhá pri štarte, ak chýba vlastnícke meno alebo heslo;
- HTTP Basic Auth rozlišuje rolu `owner` a `demo` pomocou constant-time
  porovnania;
- TLS končí na Railway edge;
- odpovede nastavujú `Cache-Control: no-store`, HSTS, `nosniff`, zákaz iframe,
  nulový referrer a zákaz kamery/mikrofónu/geolokácie;
- upload sa streamuje po 1 MB a pri prekročení limitu sa čiastočný súbor zmaže;
- názov uploadu dostane náhodný prefix a download povoľuje iba basename v
  pracovnom adresári;
- verejný režim vypína online anotáciu a BAM adresár;
- jeden procesový `JOB_LOCK` povoľuje najviac jednu imputáciu naraz;
- pracovné súbory sa mažú približne po 6 hodinách a vždy sa stratia pri
  redeployi; referenčný volume sa nemaže;
- aplikačné logy nie sú vo verejnom režime dostupné cez API.

HTTP Basic Auth je prevádzkové minimum, nie kompletný systém používateľských
účtov. Pred spracovaním genetických údajov tretích osôb treba doplniť právny
základ, informovaný súhlas, audit, incidentný proces a podľa požiadaviek aj
šifrovanie uložených dát vlastným kľúčom.

## 12. Prevádzka a diagnostika

Základná kontrola:

```bash
curl https://genome-normalizer-production.up.railway.app/healthz
```

Po vlastníckom prihlásení `/api/status` musí hlásiť `java: true`,
`beagle_jar: true`, `fasta: true`, `bref3_files: 23`, `default_panel: true` a
`imputation_ready: true`.

Lokálne testy:

```bash
python -m unittest discover -s tests -v
bash -n targeted_impute.sh
bash -n setup_reference.sh
bash -n deploy/start.sh
```

Overenie existujúcej referencie bez siete:

```bash
REFERENCE_DIR=/app/ref bash setup_reference.sh --verify-only
```

Najčastejšie poruchy:

- `imputation_ready=false`: skontrolovať Java, cesty premenných a obsah volume;
- chýbajúci chromozóm: porovnať súbor s `reference_manifest.tsv` a znovu
  spustiť bootstrap;
- HTTP 401: nesprávne Basic Auth údaje;
- HTTP 403 na DEMO účte: požadovaná funkcia je úmyselne povolená iba vlastníkovi;
- HTTP 409: iná imputácia drží globálny zámok;
- HTTP 413: upload alebo celá požiadavka prekročila nastavený limit;
- nedostatočný sibling prienik: výsledok označiť ako nedostatočné dáta, nie
  dopĺňať 561 panelových SNP do IBD vetvy.

## 13. Limity

- Základná zostava je GRCh37/hg19; hg38 vyžaduje validovaný liftover.
- Produkčné WGS/WES pomocné referencie nie sú v aktuálnom Railway nasadení
  kompletné; čipová normalizácia a cielená panelová imputácia sú pripravené.
- Presná uzavretá markerová šablóna platformy v5 nie je súčasťou repozitára.
  Predvolený výstup preto zachováva plný meraný markerový set v otvorenej
  normalizovanej schéme.
- Predispozičný výstup je výskumný a informatívny, nie klinická diagnóza.
  Klinicky významný nález musí potvrdiť akreditované laboratórium.

