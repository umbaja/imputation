FROM python:3.11-slim

# systemove zavislosti pre cielenu imputaciu:
#   default-jre-headless -> Beagle (java)
#   bcftools, tabix       -> bcftools + bgzip + tabix (htslib)
#   samtools              -> faidx index FASTA (pre bcftools +fixref)
#   curl, unzip           -> stahovanie referencie a plink2
RUN apt-get update && apt-get install -y --no-install-recommends \
    default-jre-headless bcftools tabix samtools curl unzip ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# plink2 (nie je v apt) — stiahni oficialny linux binary
RUN curl -fsSL -o /tmp/plink2.zip \
      "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip" \
    && unzip -o /tmp/plink2.zip -d /usr/local/bin \
    && chmod +x /usr/local/bin/plink2 \
    && rm /tmp/plink2.zip \
    && plink2 --version

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# referenciu (panel/beagle/fasta) pripoj ako volume do /app/ref a nastav premenne:
ENV BEAGLE_JAR=/app/ref/beagle.jar \
    REF_DIR=/app/ref/b37.bref3 \
    MAP_DIR=/app/ref/maps \
    FASTA=/app/ref/human_g1k_v37.fasta

EXPOSE 8000
CMD ["python3","-m","uvicorn","app.main:app","--host","0.0.0.0","--port","8000"]
