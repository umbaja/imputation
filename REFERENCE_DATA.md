# Referenčné dáta

Veľké referenčné súbory nie sú súčasťou Git repozitára ani Docker obrazu.
Ich presná verzia, zdroj, veľkosť a SHA-256 sú zaznamenané v
[`reference_manifest.tsv`](reference_manifest.tsv). Skript
[`setup_reference.sh`](setup_reference.sh) podľa manifestu súbory stiahne,
overí a až potom ich sprístupní aplikácii v lokálnom priečinku `ref/`.

## Obsah referencie

| Komponent | Verzia / zostava | Zdroj |
|---|---|---|
| Beagle | 5.5, `27Feb25.75f` | [Beagle, University of Washington](https://faculty.washington.edu/browning/beagle/beagle.html) |
| Imputačný panel | 1000 Genomes Phase 3 v5a, `bref3` | distribúcia Beagle, GRCh37/b37 |
| Genetické mapy | PLINK GRCh37 | distribúcia Beagle |
| Referenčná FASTA | `human_g1k_v37.fasta` | [1000 Genomes / ENA FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/) |

Pred ďalšou redistribúciou externých dát treba overiť podmienky príslušného
datasetu. Tento repozitár distribuuje iba manifest a sťahovací mechanizmus;
nevytvára kópiu panelu 1000 Genomes ani FASTA. Beagle má vlastnú licenciu
uvedenú na stránke autora.

## Inštalácia a kontrola

Kompletná referencia pre chromozómy 1–22 a X:

```bash
bash setup_reference.sh
```

Iba vybrané chromozómy, napríklad pre skúšobný beh:

```bash
CHROMS="20 21 22" bash setup_reference.sh
```

Kontrola už stiahnutých súborov bez sieťovej komunikácie:

```bash
bash setup_reference.sh --verify-only
```

Prerušené sťahovanie zostane v súbore s príponou `.part` a ďalší beh sa ho
pokúsi obnoviť. Existujúci súbor s nesprávnym checksumom sa neprepíše skôr,
než je úspešne stiahnutá a overená náhrada; následne zostane zachovaný s
príponou `.invalid-<čas>`.

Možno použiť aj iné umiestnenie dát:

```bash
REFERENCE_DIR=/mnt/genome-reference bash setup_reference.sh
```

## Aktualizácia verzie

Zmena URL alebo verzie referencie musí byť samostatná vedomá zmena manifestu.
Pri aktualizácii treba upraviť zároveň URL, počet bajtov a SHA-256 a znovu
otestovať normalizáciu aj imputáciu. Nepoužíva sa pohyblivý odkaz na „najnovšiu“
verziu, aby sa rovnaký commit dal neskôr reprodukovať.
