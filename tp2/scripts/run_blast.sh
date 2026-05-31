#!/usr/bin/env bash
set -euo pipefail

mkdir -p resultados

OUTFMT_BASE="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq stitle"
THREADS="${THREADS:-4}"
MAX_TARGET_SEQS="${MAX_TARGET_SEQS:-10}"

WORKDIR="${TMPDIR:-/tmp}/tp2_blast_work"
mkdir -p "$WORKDIR/data" "$WORKDIR/resultados"
cp data/drrA.faa data/map2k7.faa data/ndk.faa data/VFDB_setB_pro.fas data/staph_aureus_N315_proteome.faa "$WORKDIR/data/"
if [ -s data/enterotoxinas.faa ]; then
  cp data/enterotoxinas.faa "$WORKDIR/data/"
fi

makeblastdb -dbtype prot -in "$WORKDIR/data/VFDB_setB_pro.fas" -out "$WORKDIR/data/VFDB_setB_pro"

for query in drrA map2k7 ndk; do
  blastp \
    -query "$WORKDIR/data/${query}.faa" \
    -db "$WORKDIR/data/VFDB_setB_pro" \
    -evalue 1e-5 \
    -num_threads "$THREADS" \
    -max_target_seqs "$MAX_TARGET_SEQS" \
    -outfmt "$OUTFMT_BASE" \
    -out "$WORKDIR/resultados/${query}_vs_vfdb.tsv"
done

blastp \
  -query "$WORKDIR/data/staph_aureus_N315_proteome.faa" \
  -db "$WORKDIR/data/VFDB_setB_pro" \
  -evalue 1e-5 \
  -num_threads "$THREADS" \
  -max_target_seqs "$MAX_TARGET_SEQS" \
  -outfmt "$OUTFMT_BASE" \
  -out "$WORKDIR/resultados/N315_vs_VFDB.tsv"

if [ -s "$WORKDIR/data/enterotoxinas.faa" ]; then
  makeblastdb -dbtype prot -in "$WORKDIR/data/enterotoxinas.faa" -out "$WORKDIR/data/enterotoxinas"
  blastp \
    -query "$WORKDIR/data/staph_aureus_N315_proteome.faa" \
    -db "$WORKDIR/data/enterotoxinas" \
    -evalue 1e-5 \
    -num_threads "$THREADS" \
    -max_target_seqs "$MAX_TARGET_SEQS" \
    -outfmt "$OUTFMT_BASE" \
    -out "$WORKDIR/resultados/N315_vs_enterotoxinas.tsv"
fi

cp "$WORKDIR"/resultados/*.tsv resultados/
