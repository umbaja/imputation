# Bezpečnosť a genetické dáta

Genotypové súbory sú citlivé osobné a zdravotné údaje.

- Surové vzorky, VCF/BAM/PLINK súbory, logy a výsledky necommitujte do Gitu.
- Webovú aplikáciu štandardne spúšťajte iba na `127.0.0.1`.
- Priečinky `samples/`, `uploads/`, `work/`, `out/`, `results/` a `ref/` sú
  zámerne ignorované.
- Pred zdieľaním logu skontrolujte názvy súborov, identifikátory vzoriek a
  absolútne cesty.
- Verejné nasadenie vyžaduje autentifikáciu, TLS, šifrované úložisko, audit,
  časovo obmedzenú retenciu a právny základ na spracovanie genetických údajov.
- Výsledky predispozícií nie sú klinická diagnóza.
- Zdieľaný DEMO účet spracúva iba verzovaný syntetický súbor. Middleware mu
  blokuje vlastné uploady, všeobecné downloady, logy a imputáciu.
- Zdieľaný TEST účet smie používať iba kompletnú produkčnú dávkovú pipeline.
  Server vynucuje spoločný limit päť prijatých genotypov denne a počítadlo
  ukladá na persistentný Railway volume.

Bezpečnostnú chybu nezverejňujte spolu s reálnymi genotypovými dátami. Pošlite
iba minimálny syntetický reprodukčný príklad.
