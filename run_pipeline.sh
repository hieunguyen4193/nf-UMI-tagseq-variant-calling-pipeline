# Run N E X T F L O W pipeline for ctDNA project.
# ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
# Define input parameters.
# ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
RUN=H008;
mkdir -p /media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/ECD/output/$RUN;
maindir=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/ECD;
srcdir=$maindir/src_$RUN;
for group in liver;do \
CPU_BWA=5;
re=-resume;
RESOURCES=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/data/annotation_resources;
BWA_REPO=${RESOURCES}/HG38DIR/hg38_selected;
# --------
INPUT_DIR=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/data/ctDNA/fastq/R778/$group;
OUTPUT_DIR=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/ECD/output/$RUN/$group;
mkdir -p $OUTPUT_DIR;

BEDFILE_colon=${RESOURCES}/target_genes/new_version/oncoTTYSH.targets.hg38.bed;
BEDFILE_liver=${RESOURCES}/target_genes/new_version/onco08.targets.hg38.bed;
if [ "$group" = "liver" ]; then BEDFILE=$BEDFILE_liver;
BAITFILE=${RESOURCES}/target_genes/new_version/onco08.baits.hg38.interval_list;
TARGET_INTERVAL_FILE=${RESOURCES}/target_genes/new_version/onco08.targets.hg38.interval_list;
GENE_PANEL="onco08"
else BEDFILE=$BEDFILE_colon;
BAITFILE=${RESOURCES}/target_genes/new_version/oncoTTYSH.baits.hg38.interval_list;
TARGET_INTERVAL_FILE=${RESOURCES}/target_genes/new_version/oncoTTYSH.targets.hg38.interval_list;
GENE_PANEL="oncoTTYSH"
fi;
echo $BEDFILE >> $srcdir/log.txt;
echo $BAITFILE >> $srcdir/log.txt;
echo $TARGET_INTERVAL_FILE >> $srcdir/log.txt;
echo $GENE_PANEL >> $srcdir/log.txt;

# --------
BASH_DIR=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/ECD/src_H001/bash_script;

HG38DIR=$BWA_REPO.fa;


dir_cache=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/data/dir_cache_VEP99;

fasta=$dir_cache/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz;
clinvar=$dir_cache/clinvar/clinvar_20201026.vcf.gz;
clinvar_index=$dir_cache/clinvar/clinvar_20201026.vcf.gz.tbi
COSMIC=$dir_cache/COSMIC/CosmicCodingMuts.normal.vcf.gz;
cosmic_index=$dir_cache/COSMIC/CosmicCodingMuts.normal.vcf.gz.tbi;

nextflow=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/data/nextflow_exe/nextflow;

chmod -R a+rwx $RESOURCES;
chmod -R a+rwx $dir_cache;
chmod -R a+rwx $OUTPUT_DIR;

$nextflow $srcdir/ECD_pipeline_tagseq.nf -c $srcdir/ECD_pipeline_tagseq.config \
-w $srcdir/work \
--BWA_REPO $BWA_REPO \
--CPU_BWA $CPU_BWA \
--INPUT_DIR $INPUT_DIR \
--OUTPUT_DIR $OUTPUT_DIR \
--BASH_DIR $BASH_DIR \
--BEDFILE $BEDFILE \
--HG38DIR $HG38DIR \
--RESOURCES $RESOURCES \
--BAIT $BAITFILE \
--GENE_PANEL $GENE_PANEL \
--TARGET $TARGET_INTERVAL_FILE \
--dir_cache $dir_cache \
--fasta $fasta \
--clinvar_index $clinvar_index \
--clinvar $clinvar \
--COSMIC $COSMIC \
--COSMIC_INDEX $cosmic_index \
$re;done