# Día 1
## Uso de tablas para el análisis de datos en serie
| Sample_ID       | raw_read_R1 | raw_read_R2 | clean_read_R1 | clean_read_R2 | assembly |
|-----------------|-------------|-------------|---------------|---------------|----------|
| P3_2m_110220_MG | /home/wcastillo/01_raw_reads/P3_2m_110220_gi_CTTAAGTGAC-TACGCTAGTT_L001_R1_001.fastq.gz | /home/wcastillo/01_raw_reads/P3_2m_110220_gi_CTTAAGTGAC-TACGCTAGTT_L001_R2_001.fastq.gz | /home/wcastillo/02_clean_reads/P3_2m_110220_MG_clean_R1.fastq.gz | /home/wcastillo/02_clean_reads/P3_2m_110220_MG_clean_R2.fastq.gz | /home/wcastillo/03_assembly/P3_2m_110220_MG/final.contigs.fa |
| P3_2m_210122_MG | /home/wcastillo/01_raw_reads/P3_2m_210122_gi_CGTGGAACAC-GAACAGATGG_L002_R1_001.fastq.gz | /home/wcastillo/01_raw_reads/P3_2m_210122_gi_CGTGGAACAC-GAACAGATGG_L002_R2_001.fastq.gz | /home/wcastillo/02_clean_reads/P3_2m_210122_MG_clean_R1.fastq.gz | /home/wcastillo/02_clean_reads/P3_2m_210122_MG_clean_R2.fastq.gz | /home/wcastillo/03_assembly/P3_2m_210122_MG/final.contigs.fa |
| P3_2m_280122_MG | /home/wcastillo/01_raw_reads/P3_2m_280122_gi_ATGGAAGTGG-CACCTACGAA_L002_R1_001.fastq.gz | /home/wcastillo/01_raw_reads/P3_2m_280122_gi_ATGGAAGTGG-CACCTACGAA_L002_R2_001.fastq.gz | /home/wcastillo/02_clean_reads/P3_2m_280122_MG_clean_R1.fastq.gz | /home/wcastillo/02_clean_reads/P3_2m_280122_MG_clean_R2.fastq.gz | /home/wcastillo/03_assembly/P3_2m_280122_MG/final.contigs.fa |
| P3_2m_140122_MG | /home/wcastillo/01_raw_reads/P3_2m_140122_gi_ACTAGTGCTT-ACTCTGTTCT_L002_R1_001.fastq.gz | /home/wcastillo/01_raw_reads/P3_2m_140122_gi_ACTAGTGCTT-ACTCTGTTCT_L002_R2_001.fastq.gz | /home/wcastillo/02_clean_reads/P3_2m_140122_MG_clean_R1.fastq.gz | /home/wcastillo/02_clean_reads/P3_2m_140122_MG_clean_R2.fastq.gz | /home/wcastillo/03_assembly/P3_2m_140122_MG/final.contigs.fa |
| P3_30m_140122_MG | /home/wcastillo/01_raw_reads/P3_30m_140122_gi_TGAGATCACA-TGTACCGTGC_L002_R1_001.fastq.gz | /home/wcastillo/01_raw_reads/P3_30m_140122_gi_TGAGATCACA-TGTACCGTGC_L002_R2_001.fastq.gz | /home/wcastillo/02_clean_reads/P3_30m_140122_MG_clean_R1.fastq.gz | /home/wcastillo/02_clean_reads/P3_30m_140122_MG_clean_R2.fastq.gz | /home/wcastillo/03_assembly/P3_30m_140122_MG/final.contigs.fa |


## ¿Para qué tabular rutas?

```
#!/bin/bash
#SBATCH --job-name=stats_assemblies
#SBATCH --output=logs_stats/stats_%j.out
#SBATCH --error=logs_stats/stats_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00

METADATA="antarctica_full_paths.tsv"
TMPDIR="./tmp_stats"
OUTPUT_STATS="assembly_read_stats_full.tsv"
FILTERED_DIR="./filtered_contigs"
COMBINED_FASTA="UViGs_combined.fasta"

mkdir -p "$TMPDIR" "$FILTERED_DIR"
> "$OUTPUT_STATS"
> "$COMBINED_FASTA"

# Write header
echo -e "Sample_ID\tRaw_Reads_Total_bp\tClean_Reads_Total_bp\tNum_Contigs\tTotal_Contigs_bp\tLongest_Contig\tN50\tNum_Contigs_>10kb\tN50_>10kb" > "$OUTPUT_STATS"

tail -n +2 "$METADATA" | while IFS=$'\t' read -r SAMPLE R1 R2 CLEAN_R1 CLEAN_R2 ASM; do
    echo "Processing $SAMPLE..."

    # ---- Read stats ----
    RAW_TOTAL=$( (zcat "$R1" "$R2" 2>/dev/null || cat "$R1" "$R2") | awk 'NR % 4 == 2 {sum += length($0)} END {print sum}')
    CLEAN_TOTAL=$( (zcat "$CLEAN_R1" "$CLEAN_R2" 2>/dev/null || cat "$CLEAN_R1" "$CLEAN_R2") | awk 'NR % 4 == 2 {sum += length($0)} END {print sum}')

    # ---- Assembly stats ----
    NUM_CONTIGS=$(grep -c "^>" "$ASM")
    CONTIG_LENGTHS=($(awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' "$ASM" | sort -nr))
    TOTAL_LEN=0
    for len in "${CONTIG_LENGTHS[@]}"; do ((TOTAL_LEN+=len)); done
    LONGEST="${CONTIG_LENGTHS[0]}"
    HALF_TOTAL=$((TOTAL_LEN / 2))
    N50=0; CUM=0
    for len in "${CONTIG_LENGTHS[@]}"; do
        ((CUM+=len))
        if ((CUM >= HALF_TOTAL)); then N50=$len; break; fi
    done

    # ---- Stats for contigs >10kb ----
    CONTIGS_OVER_10K=()
    for len in "${CONTIG_LENGTHS[@]}"; do
        if (( len > 10000 )); then CONTIGS_OVER_10K+=("$len"); fi
    done
    NUM_OVER_10K=${#CONTIGS_OVER_10K[@]}
    TOTAL_OVER_10K=0
    for len in "${CONTIGS_OVER_10K[@]}"; do ((TOTAL_OVER_10K+=len)); done
    HALF_OVER_10K=$((TOTAL_OVER_10K / 2))
    N50_OVER_10K=0; CUM=0
    for len in "${CONTIGS_OVER_10K[@]}"; do
        ((CUM+=len))
        if ((CUM >= HALF_OVER_10K)); then N50_OVER_10K=$len; break; fi
    done

    # ---- Write line ----
    echo -e "$SAMPLE\t$RAW_TOTAL\t$CLEAN_TOTAL\t$NUM_CONTIGS\t$TOTAL_LEN\t$LONGEST\t$N50\t$NUM_OVER_10K\t$N50_OVER_10K" >> "$OUTPUT_STATS"

done  # ← this was missing

echo "[DONE] Stats written to: $OUTPUT_STATS"
```
