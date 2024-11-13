#!/bin/bash
#SBATCH --job-name=d9_bedtools                          # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=esro2814@colorado.edu        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:10:00                                 # Time limit hrs:min:sec
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH -p short
#SBATCH --output=/scratch/Users/esro2814/workshop-day9/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/esro2814/workshop-day9/eofiles/%x_%j.err


################## USER INPUT ######################################

chrom_sizes='/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.chrom.sizes'
ANNOT='/scratch/Shares/public/sread2024/cookingShow/day9b/annotation/hg38_refseq_cleaned.sorted.bed'
CHIP_DIR='/scratch/Users/esro2814/workshop-day8/macs/'
HCT116_DMSO_PEAKS=${CHIP_DIR}/HCT116_DMSO_peaks_clean.narrowPeak
HCT116_NUTLIN_PEAKS=${CHIP_DIR}/HCT116_NUTLIN_peaks_clean.narrowPeak
results='/scratch/Users/esro2814/workshop-day9/bedtools_results/'
TSS_ANNOT=${results}/hg38_refseq_tssregion.bed

mkdir -p ${results}
################## JOB INFO ########################################

printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

################### LOAD MODULES ###################################

module load bedtools/2.28.0


################## SCRIPT ##########################################

# Write your bedtools/file manipulation commands here

# --- Which p53 peaks are associated with genes (not intergenic)
bedtools intersect -wa -a $HCT116_DMSO_PEAKS -b $ANNOT > ${results}p53_DMSO_peaks_in_genes.bed
bedtools intersect -wa -a $HCT116_NUTLIN_PEAKS -b $ANNOT > ${results}p53_Nutlin_peaks_in_genes.bed
# --- Which p53 peaks are entirely within genes 
bedtools intersect -wa -f 1 -a $HCT116_DMSO_PEAKS -b $ANNOT > ${results}DMSO_p53_peaks_entirely_in_genes.bed
bedtools intersect -wa -f 1 -a $HCT116_NUTLIN_PEAKS -b $ANNOT > ${results}Nutlin_p53_peaks_entirely_in_genes.bed
# --- Which p53 peaks are intergenic 
bedtools intersect -wa -v -a $HCT116_DMSO_PEAKS -b $ANNOT > ${results}intergenic_DMSO_p53_peaks.bed
bedtools intersect -wa -v -a $HCT116_NUTLIN_PEAKS -b $ANNOT > ${results}intergenic_Nutlin_p53_peaks.bed
# --- Which genes overlap with p53 peaks? 
bedtools intersect -wa -F 0.10  -a $ANNOT -b $HCT116_DMSO_PEAKS > ${results}genes_with_p53_peaks_DMSO.bed
bedtools intersect -wa -F 0.10  -a $ANNOT -b $HCT116_NUTLIN_PEAKS > ${results}genes_with_p53_peaks_Nutlin.bed

# --- Which genes have p53 ChIP sites in their promoters

cat "ANNOT" \
| awk '{FS=OFS="\t"}; {if ($6 == "+") print $1,$2,($2+1),$4,".",$6}' \
> "$results"/hg38_refseq_tss_pos.bed
cat "$ANNOT" \
| awk '{FS=OFS="\t"}; {if ($6 == "-") print $1,($3-1),$3,$4,".",$6}' \
> "$results"/hg38_refseq_tss_neg.bed
cat "$results"/hg38_refseq_tss_pos.bed \
"$results"/hg38_refseq_tss_neg.bed \
> "$results"/hg38_refseq_tss.bed

###---- Now do for Nutlin treated HCT116 Samples ----###




###---- Jacquard Calculations ----###
