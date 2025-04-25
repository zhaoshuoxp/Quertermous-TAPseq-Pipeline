
# TAP-seq + SCEPTRE Analysis Pipeline

This repository contains a streamlined pipeline to perform TAP-seq alignment with **Cell Ranger**, generate filtered feature-barcode matrices, and perform **SCEPTRE** single-cell perturbation analysis.

## 📁 Directory Structure

```
├── run_pipeline.sh       # Main bash script (build reference, align, call R analysis)
├── sceptre.R             # R script (analyze filtered matrices with SCEPTRE)
├── refdata-cellranger-custom/   # (Auto-generated) Custom reference genome for TAP-seq
├── all_sample/           # (Auto-generated) Cell Ranger alignment output
├── sceptre/              # (Auto-generated) SCEPTRE analysis output
├── [Input files]         # Provided input files
│   ├── genelist.txt
│   ├── guides.txt
│   ├── lib.csv
│   ├── grna_target.txt
│   ├── feature_large.csv
│   ├── pos_controls.txt
│   └── discovery_cis.txt
└── README.md             # This instruction file
```

## 📋 Requirements

- `cellranger8`
- `R (>= 4.1)` with packages:
  - `sceptre`
  - `Seurat`
  - `reshape2`

The bash script will automatically check for `cellranger8`.  
The R script will automatically install missing R packages if needed.

## ⚙️ Input Files

| Filename | Description |
| :--- | :--- |
| `genelist.txt` | List of target genes for GTF filtering |
| `guides.txt` | List of sgRNA guide sequences captured |
| `lib.csv` | Cell Ranger libraries CSV file |
| `grna_target.txt` | Mapping between sgRNA and target regions |
| `feature_large.csv` | Feature reference for Cell Ranger |
| `pos_controls.txt` | Positive control grna-target pairs for QC |
| `discovery_cis.txt` | Discovery grna-target pairs for analysis |

## 🚀 How to Run

1. **Prepare all input files** and place them in the working directory.
2. **Ensure** `run_pipeline.sh` and `sceptre.R` are executable:
   ```bash
   chmod +x run_pipeline.sh sceptre.R
   ```
3. **Run the full pipeline** with:
   ```bash
   ./run_pipeline.sh <fastq_dir> <ref_path> <genelist.txt> <feature_large.csv> <lib.csv> <guides.txt> <covariate_file> <grna_target.txt> <pos_controls.txt> <discovery_cis.txt>
   ```

   Example:
   ```bash
   ./run_pipeline.sh fastq/ /path/to/refdata/ genelist.txt feature_large.csv lib.csv guides.txt covariate.tsv grna_target.txt pos_controls.txt discovery_cis.txt
   ```

## 🔎 Pipeline Details

1. **Reference Generation**:
   - Filter GTF to keep only `protein_coding` genes.
   - Further filter to only retain target genes from `genelist.txt`.
   - Build custom Cell Ranger reference.

2. **Cell Ranger Count**:
   - Align reads using the custom reference genome.
   - Generate gene expression and CRISPR feature matrices.

3. **SCEPTRE Analysis**:
   - Import filtered matrices.
   - Assign gRNAs.
   - Run quality control, calibration check, and power check.
   - Perform discovery analysis.
   - Output results (`results.txt`) and save full object (`sceptre.rds`).

## 📄 Outputs

After successful execution:

| File/Folder | Description |
| :--- | :--- |
| `all_sample/` | Cell Ranger filtered output |
| `sceptre/` | SCEPTRE results and figures |
| `sceptre/results.txt` | Discovery analysis results (with adjusted p-values) |
| `sceptre/sceptre.rds` | Full R object of the analysis |

## 🛠️ Troubleshooting

- Make sure `cellranger8` is correctly loaded (e.g., `module load cellranger/8.0.0` if using HPC).
- Ensure R packages can be installed from CRAN or Bioconductor.
- Check available memory (200GB recommended) and threads (default 20 cores).
- Provide all 10 required input parameters when running `run_pipeline.sh`.

## ✨ Notes

- This pipeline is optimized for **TAP-seq CRISPR guide capture** experiments.
- You can customize memory, thread number, or paths by modifying the top of `run_pipeline.sh`.
- All outputs will be generated in the working directory by default.

# 📬 Contact

For questions or issues, please contact the pipeline maintainer.
