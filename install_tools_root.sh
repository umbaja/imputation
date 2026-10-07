#!/usr/bin/env bash
# Spusti AKO ROOT (v PowerShelli: wsl -u root, potom tento skript).
# Nainstaluje systemovo java, bcftools, samtools, tabix, plink2 pre vsetkych.
set -e
apt-get update -y
apt-get install -y default-jre bcftools tabix samtools wget curl unzip python3 python3-venv python3-pip
# plink2 systemovo do /usr/local/bin
if ! command -v plink2 >/dev/null 2>&1; then
  wget -q https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64.zip -O /tmp/plink2.zip
  unzip -o /tmp/plink2.zip -d /usr/local/bin/
  chmod +x /usr/local/bin/plink2
fi
# plink 1.9 (ma --23file, potrebne na prevod 23andMe -> VCF)
if ! command -v plink >/dev/null 2>&1; then
  wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip -O /tmp/plink1.zip
  unzip -o /tmp/plink1.zip plink -d /usr/local/bin/
  chmod +x /usr/local/bin/plink
fi
echo "=== hotovo, verzie: ==="
java -version; bcftools --version | head -1; plink2 --version; samtools --version | head -1
echo ">> Nastroje nainstalovane. Vrat sa do bezneho pouzivatelskeho terminalu a spusti setup_reference.sh."
