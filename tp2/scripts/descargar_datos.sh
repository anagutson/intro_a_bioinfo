#!/usr/bin/env bash
set -euo pipefail

mkdir -p data resultados

curl -L "https://rest.uniprot.org/uniprotkb/Q29ST3.fasta" -o data/drrA.faa
curl -L "https://www.uniprot.org/uniprot/O14733.fasta" -o data/map2k7.faa
curl -L "https://www.uniprot.org/uniprot/P68870.fasta" -o data/ndk.faa

# En algunas redes el nombre www.mgc.ac.cn no resuelve, pero el servidor VFDB
# responde usando su IP con el Host header oficial.
curl -L --header "Host: www.mgc.ac.cn" "http://183.242.79.135/VFs/Down/VFDB_setB_pro.fas.gz" -o data/VFDB_setB_pro.fas.gz
gunzip -fk data/VFDB_setB_pro.fas.gz

curl -L "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?id=NC_002745.2&db=nuccore&report=fasta_cds_aa&retmode=text" -o data/staph_aureus_N315_proteome.faa

curl -L "https://docs.google.com/document/d/11mIjHAc9AWGzf8cPbISiRxpfd_KOZqb_Hff2_oJivuU/export?format=txt" -o data/enterotoxinas_raw.txt
