#!/bin/bash
set -euo pipefail

# -------------------------
required_tools=("cellranger")
for tool in "${required_tools[@]}"; do
  if ! command -v "$tool" &>/dev/null; then
    echo "[ERROR] Required tool '$tool' not found in PATH. Please install or load it first."
    exit 1
  fi
done

# ----------------------
MEMORY=200
THREADS=20
R_SCRIPT_PATH="./sceptre.R"  

# --------------------------
if [ "$#" -lt 10 ]; then
  echo "[USAGE] $0 <fastq_dir> <ref_path> <target_gene_list> <feature_ref_csv> <library_csv> <guides_file> <covariate_file> <grna_target_file> <positive_control_pairs_file> <discovery_pairs_file>"
  exit 1
fi

FASTQ_DIR="$1"
REF_PATH="$2"
TARGET_GENE_LIST="$3"
FEATURE_REF="$4"
LIBRARY_CSV="$5"
GUIDES_FILE="$6"
COV_FILE="$7"
GRNA_TARGET_FILE="$8"
POSITIVE_CONTROL_PAIRS="$9"
DISCOVERY_PAIRS="${10}"

# --------------------
echo "[INFO] Building custom Cell Ranger reference..."
REF_DIR="refdata-cellranger-custom"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

ln -sf "${REF_PATH}/genes/genes.gtf.gz" .
ln -sf "${REF_PATH}/fasta/genome.fa" .

gunzip -c genes.gtf.gz > genes.gtf

cellranger mkgtf genes.gtf genes.filter.gtf --attribute=gene_type:protein_coding

grep -w -f "$TARGET_GENE_LIST" genes.filter.gtf > genes.filter2.gtf

cellranger mkref --genome=TAPseq_genome --fasta=genome.fa --genes=genes.filter2.gtf
cd ..

# -----------------------
echo "[INFO] Running Cell Ranger count..."
cellranger count \
  --localcores="$THREADS" \
  --localmem="$MEMORY" \
  --id="all_sample" \
  --create-bam=true \
  --libraries="$LIBRARY_CSV" \
  --transcriptome="${REF_DIR}/TAPseq_genome" \
  --feature-ref="$FEATURE_REF"

# -------------------
echo "[INFO] Running SCEPTRE analysis with R script..."

MTX_DIR="all_sample/outs/filtered_feature_bc_matrix"

mkdir -p sceptre
cd sceptre

chmod +x "$R_SCRIPT_PATH"
"$R_SCRIPT_PATH" "../${MTX_DIR}" "$GUIDES_FILE" "$COV_FILE" "$GRNA_TARGET_FILE" "$POSITIVE_CONTROL_PAIRS" "$DISCOVERY_PAIRS" "$THREADS"

cd ..

echo "[INFO] All steps completed successfully!"
