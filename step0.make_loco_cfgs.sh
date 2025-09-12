#!/usr/bin/env bash
set -euo pipefail

GWAS_DIR="/data/wangdi/CGRProjects/sonja/DEPICT/run30gwas/data/gwas"
CFG_DIR="/data/wangdi/CGRProjects/sonja/DEPICT/run30gwas/configs_loco"
TEMPLATE="/data/wangdi/CGRProjects/sonja/DEPICT/run30gwas/depict.loco.template.cfg"

mkdir -p "$CFG_DIR"

> "$CFG_DIR/configs.list2"

for f in "$GWAS_DIR"/*.txt; do
  gwaslabel=$(basename "$f" .txt)
  for chr in {1..22}; do
    outcfg="$CFG_DIR/${gwaslabel}.chr${chr}.cfg"
    sed -e "s|{GWASFILE}|$f|g" \
        -e "s|{GWASLABEL}|$gwaslabel|g" \
        -e "s|{CHR}|$chr|g" \
        "$TEMPLATE" > "$outcfg"
    echo "$outcfg" >> "$CFG_DIR/configs.list2"
  done
done

echo "Wrote $(wc -l < "$CFG_DIR/configs.list2") configs."

