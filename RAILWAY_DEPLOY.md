# Nasadenie mimo lokálneho počítača

Genome Normalizer je navrhnutý predovšetkým ako lokálna aplikácia. Genotypové
dáta sú citlivé a plná imputácia potrebuje približne 13 GB referencií.

## Odporúčanie

Pre jedného používateľa používajte WSL/Linux alebo Docker na lokálnom počítači.
Pre tím použite interný server s autentifikáciou, TLS, šifrovaným diskom a
časovo obmedzeným mazaním uploadov.

## Railway a podobné PaaS

Docker image možno zostaviť z koreňa repozitára, ale referencie nesmú byť jeho
súčasťou. Služba potrebuje persistentný volume pripojený ako `/app/ref` a tieto
premenné:

```text
BEAGLE_JAR=/app/ref/beagle.jar
REF_DIR=/app/ref/b37.bref3
MAP_DIR=/app/ref/maps
FASTA=/app/ref/human_g1k_v37.fasta
```

Pred verejným nasadením treba doplniť autentifikáciu a pravidlá retencie.
Súčasné FastAPI rozhranie ich neposkytuje, preto ho nevystavujte priamo na
verejný internet. Samotný GitHub repozitár neobsahuje genetické vzorky ani
veľké referenčné dáta.
