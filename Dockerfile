FROM python:3.11-slim

# systemove zavislosti pre cielenu imputaciu:
#   default-jre-headless -> Beagle (java)
#   bcftools, tabix       -> bcftools + bgzip + tabix (htslib)
#   samtools              -> faidx index FASTA (pre bcftools +fixref)
#   curl, unzip, gzip     -> overené stiahnutie referencií pri prvom štarte
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash coreutils gawk gzip tzdata default-jre-headless bcftools tabix samtools \
    curl unzip ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# referenciu (panel/beagle/fasta) pripoj ako volume do /app/ref a nastav premenne:
ENV BEAGLE_JAR=/app/ref/beagle.jar \
    REF_DIR=/app/ref/b37.bref3 \
    MAP_DIR=/app/ref/maps \
    FASTA=/app/ref/human_g1k_v37.fasta \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

EXPOSE 8000
CMD ["bash", "deploy/start.sh"]
