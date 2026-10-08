# Railway nasadenie

Toto nasadenie je určené pre jednu heslom chránenú inštanciu s vlastníckym,
voliteľným obmedzeným DEMO a kvótovaným TEST účtom. Genetické
vstupy sa spracujú v kontajneri, výsledky sa dajú stiahnuť a pracovné súbory
sa automaticky mažú. Imputácia beží vždy iba pre jednu úlohu naraz.

## Potrebné zdroje

Plná referencia má po stiahnutí a rozbalení približne 12–13 GB. Railway
služba preto potrebuje persistentný volume s rezervou, odporúčane **15 GB**,
pripojený na `/app/ref`. Takáto kapacita môže vyžadovať platený Railway plán;
aktuálne limity treba overiť v [Railway pricing](https://railway.com/pricing).

Odporúčaná minimálna pamäť pre jednu imputáciu je 4 GB. Premenné `JOBS=1` a
`XMX_JOB=3g` zámerne obmedzujú paralelizmus a Java heap.

## Postup v Railway

1. Vytvorte **New Project → Deploy from GitHub repo** a vyberte
   `umbaja/imputation`.
2. K službe pridajte persistentný volume s mount path `/app/ref` a kapacitou
   aspoň 15 GB.
3. V časti **Variables** nastavte:

   ```text
   PUBLIC_MODE=1
   APP_USERNAME=<zvolené používateľské meno>
   APP_PASSWORD=<dlhé náhodné heslo>
   DEMO_USERNAME=demo
   DEMO_PASSWORD=genome-demo-2026
   TEST_USERNAME=test
   TEST_PASSWORD=<dlhé náhodné heslo odovzdané iba testerom>
   TEST_DAILY_LIMIT=5
   TEST_QUOTA_TIMEZONE=Europe/Bratislava
   TEST_QUOTA_FILE=/app/ref/test_account_quota.json
   BOOTSTRAP_REFERENCE=1
   REFERENCE_DIR=/app/ref
   BEAGLE_JAR=/app/ref/beagle.jar
   REF_DIR=/app/ref/b37.bref3
   MAP_DIR=/app/ref/maps
   FASTA=/app/ref/human_g1k_v37.fasta
   JOBS=1
   XMX_JOB=3g
   RUNTIME_RETENTION_HOURS=6
   MAX_UPLOAD_MB=100
   MAX_REQUEST_MB=250
   MAX_BATCH_SAMPLES=5
   ```

4. V **Settings → Healthcheck** nastavte cestu `/healthz` a timeout aspoň
   3600 sekúnd (alebo pridajte `RAILWAY_HEALTHCHECK_TIMEOUT_SEC=3600`).
5. Spustite redeploy. Pri prvom štarte sa referencia stiahne na volume a každý
   súbor sa overí podľa veľkosti a SHA-256 z `reference_manifest.tsv`. Prvý
   štart môže trvať desiatky minút; ďalšie štarty používajú volume.
6. Po úspešnom health checku `/healthz` zvoľte **Settings → Networking →
   Generate Domain**.
7. Otvorte doménu a prihláste sa hodnotami `APP_USERNAME` a `APP_PASSWORD`.

Zdieľané konto `DEMO_USERNAME` môže spustiť iba zabudovanú syntetickú ukážku.
Middleware mu blokuje uploady, všeobecné downloady, logy aj výpočtovo nákladnú
imputáciu. Ak DEMO účet nechcete, obe premenné `DEMO_*` vynechajte.

Zdieľané konto `TEST_USERNAME` môže spustiť ostrú kompletnú dávkovú pipeline.
Server počíta prijaté vzorky a povoľuje najviac `TEST_DAILY_LIMIT` genotypov za
deň, spoločne pre celý účet. Počítadlo sa obnovuje o polnoci v
`TEST_QUOTA_TIMEZONE` a súbor `TEST_QUOTA_FILE` zostáva na persistentnom
volume. Ak TEST účet nechcete, vynechajte `TEST_USERNAME` a `TEST_PASSWORD`;
limit, časové pásmo a cestu možno ponechať.

`APP_PASSWORD` ani `TEST_PASSWORD` nevkladajte do Git repozitára, dokumentácie
alebo screenshotov. Railway Variables sú jediné miesto, kde majú byť uložené.

## Prevádzkové správanie

- `/healthz` je verejný len pre Railway health check; ostatné cesty vyžadujú
  HTTP Basic prihlásenie.
- OpenAPI/Swagger a aplikačné logy nie sú vo verejnom režime dostupné.
- Maximálna veľkosť jedného uploadu je 100 MB a dávka má najviac 5 vzoriek.
- Obsah `/app/work` sa automaticky odstraňuje po 6 hodinách. Nie je na volume,
  takže zmizne aj pri redeployi/reštarte kontajnera.
- V sibling výstupe ostávajú iba priamo namerané genotypy. Imputované hodnoty
  sú iba v samostatnom predispozičnom výstupe a panelovom CSV.
- Verejný režim vypína online anotáciu aj vstup z BAM priečinka.

## Bez imputácie

Na lacný test používateľského rozhrania možno dočasne nastaviť
`BOOTSTRAP_REFERENCE=0` a nepripájať veľký volume. Konverzia formátov bude
fungovať, ale tlačidlo celej pipeline správne oznámi, že imputácia nie je
pripravená.
