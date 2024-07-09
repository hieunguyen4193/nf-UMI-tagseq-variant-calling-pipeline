// -----------------------------------------------------------------
// N E X T F L O W pipeline for data with tag-seq UMI.
// -----------------------------------------------------------------

params.CPU_BWA = ""
params.BWA_REPO= "" 
params.INPUT_DIR= ""
params.OUTPUT_DIR= ""
params.BEDFILE=""
bedfile=file(params.BEDFILE)
params.HG38DIR=""
hg38=file(params.HG38DIR)

params.RESOURCES=""
params.dbsnp="${params.RESOURCES}/genomes/dbSNP_INDELS/dbsnp_146.hg38.vcf"
params.ADAPTER="${params.RESOURCES}/adapters/TruSeq3-PE-2.fa"
ADAPTER=file(params.ADAPTER)
dbsnp=file(params.dbsnp)
params.AF_THR="0.00001"
params.INPUT_PAIRS= "$params.INPUT_DIR/*_R{1,2}*" 
params.dir_cache=""
dir_cache=file(params.dir_cache)
params.fasta=""
fasta=file(params.fasta)
params.clinvar=""
clinvar=file(params.clinvar)
params.clinvar_index=""
clinvar_index=file(params.clinvar_index)
params.COSMIC=""
COSMIC=file(params.COSMIC)
params.COSMIC_INDEX=""
COSMIC_INDEX=file(params.COSMIC_INDEX)
params.GENE_PANEL=""
params.BAIT=""
params.TARGET=""

params.bwa2="/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/ECD/src_H008/bwa-mem2-2.0pre2_x64-linux"
bait_intervals=file(params.BAIT)
target_intervals=file(params.TARGET)
// ----- ----- ----- CHANNEL ----- ----- -----
Channel
    .fromFilePairs( params.INPUT_PAIRS )
    .ifEmpty { error "Cannot find any reads matching: ${params.INPUT_PAIRS}"  }
    // .view()
    .into {fastqc_ch; trim_ch}

// -------------------------------------------
// PREPROCESSING UMIs
// -------------------------------------------

process fastqc_preprocess {
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/0_fastqc_results", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    // maxForks 5

    input:
        tuple sample_id, file(READS) from fastqc_ch 
    output:
        tuple sample_id, '*_fastqc.{zip,html}' into fastqc_results1
    script:
    """
    fastqc --quiet --outdir . ${READS[0]} ${READS[1]}
    """
}
process adapter_trimming {
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/1_trimmed_fastq", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1

    input:
        tuple sample_id, file(READS) from trim_ch
    output:
        tuple sample_id, "${sample_id}_R1_paired.fastq.gz", "${sample_id}_R2_paired.fastq.gz" into fastq2unmappedBAM

    script:
    """
    trimmomatic PE -phred33 -threads 25 ${READS[0]} ${READS[1]} \
    ${sample_id}_R1_paired.fastq.gz ${sample_id}_R1_unpaired.fastq.gz \
    ${sample_id}_R2_paired.fastq.gz ${sample_id}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTER}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 CROP:75 MINLEN:36
    """
}

process fastq_to_unmappedBAM{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/2_unmapped_bam", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1

    input:
        tuple sample_id, "${sample_id}_R1_paired.fastq.gz", "${sample_id}_R2_paired.fastq.gz" from fastq2unmappedBAM
    output:
        tuple sample_id, "unmapped_${sample_id}.bam" into bam2fastq, unmappedBAM
    script:
    """
    fgbio FastqToBam -i ${sample_id}_R1_paired.fastq.gz ${sample_id}_R2_paired.fastq.gz \
    -r 6M1S+T 6M1S+T -o unmapped_${sample_id}.bam -s true \
    --sample ${sample_id} --library "tag-seq" --read-group-id ${sample_id};
    """
}

process bam_to_noUMIfastq{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/3_noUMI_fastq", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1

    input: 
        tuple sample_id, "unmapped_${sample_id}.bam" from bam2fastq
    output:
        tuple sample_id, "${sample_id}_R1.fastq.gz", "${sample_id}_R2.fastq.gz" into noUMIfastq, bwamem2test
    script:
    """
    picard SamToFastq INPUT=unmapped_${sample_id}.bam \
    FASTQ=${sample_id}_R1.fastq SECOND_END_FASTQ=${sample_id}_R2.fastq;
    bgzip -c ${sample_id}_R1.fastq > ${sample_id}_R1.fastq.gz;
    bgzip -c ${sample_id}_R2.fastq > ${sample_id}_R2.fastq.gz; 
    """
}

process align_noUMI_FASTQ{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/4_align_noUMI", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1

    input:
        tuple sample_id, "${sample_id}_R1.fastq.gz", "${sample_id}_R2.fastq.gz" from noUMIfastq
    output:
        tuple sample_id, "noUMI_${sample_id}.sam" into sortBAM
    script:
    """
    bwa mem -t ${params.CPU_BWA} -M -R \
            '@RG\\tID:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA\\tPM:MINISEQ\\tSM:${sample_id}' \\
            ${params.BWA_REPO}\
            ${sample_id}_R1.fastq.gz ${sample_id}_R2.fastq.gz > noUMI_${sample_id}.sam

    """
}
process sort_convert2BAM{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/5_sort_BAM", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    
    input:
        tuple sample_id, "noUMI_${sample_id}.sam" from sortBAM
    output:
        tuple sample_id, "sorted_${sample_id}.bam" into mergeUMI
    script:
    """
    picard SortSam INPUT=noUMI_${sample_id}.sam OUTPUT=sorted_${sample_id}.sam SORT_ORDER=queryname;
    samtools view -S -b sorted_${sample_id}.sam > sorted_${sample_id}.bam;
    """
}

process mergeUMI{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/6_mergedBAM_UMI", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1

    input:
        tuple sample_id, "sorted_${sample_id}.bam", "unmapped_${sample_id}.bam" from mergeUMI.join(unmappedBAM)
    output:
        tuple sample_id, "mapped_withUMI_${sample_id}.bam", "mapped_withUMI_${sample_id}.bai"  into group_reads_UMI
    script:
    """
    picard MergeBamAlignment -ALIGNED sorted_${sample_id}.bam \
    -UNMAPPED unmapped_${sample_id}.bam \
    -OUTPUT mapped_withUMI_${sample_id}.bam \
    -REFERENCE_SEQUENCE ${params.BWA_REPO}.fa \
    -SORT_ORDER coordinate \
    -ALIGNER_PROPER_PAIR_FLAGS true \
    -ALIGNED_READS_ONLY true \
    -CREATE_INDEX true \
    -VALIDATION_STRINGENCY SILENT \
    -EXPECTED_ORIENTATIONS FR \
    -MAX_INSERTIONS_OR_DELETIONS -1;
    
    """
}

process group_reads_by_UMI{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/7_group_reads_byUMI", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1

    input:
        tuple sample_id, "mapped_withUMI_${sample_id}.bam", "mapped_withUMI_${sample_id}.bai" from group_reads_UMI 
    output:
        tuple sample_id, "groupedUMI_${sample_id}.bam", "${sample_id}_familySize.txt" into callconsensus
    script:
    """
    fgbio GroupReadsByUmi -i mapped_withUMI_${sample_id}.bam -o groupedUMI_${sample_id}.bam -s paired -f ${sample_id}_familySize.txt;
    """
}

process CallMolecularConsensus {
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/8_consensusBAM", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1

    input:
        tuple sample_id, "groupedUMI_${sample_id}.bam" from callconsensus
    output:
        tuple sample_id, "consensus_${sample_id}.sorted.bam" into generate_fastq
    script:
    """
    fgbio CallMolecularConsensusReads -M 3 \
    -i groupedUMI_${sample_id}.bam \
    -o ${sample_id}_consensus.bam;
    samtools sort -n ${sample_id}_consensus.bam -o consensus_${sample_id}.sorted.bam
    """
}
// -----------------------------------------
// POST-UMI PROCESSING - CALLING VARIANTS
// -----------------------------------------
process Consensus_BAM_to_FASTQ{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/9_consensusFASTQ", mode: "copy"
    // errorStrategy 'retry'
    // maxRetries 1
    input:
        tuple sample_id, "consensus_${sample_id}.sorted.bam" from generate_fastq
    output:
        tuple sample_id, "consensus_${sample_id}_R1.fastq.gz", "consensus_${sample_id}_R2.fastq.gz" into fastq_postconsensus, fastq_QC2, bwamem2_test2
    script:
    """
    bedtools bamtofastq -i consensus_${sample_id}.sorted.bam -fq consensus_${sample_id}_R1.fastq -fq2 consensus_${sample_id}_R2.fastq;
    bgzip -c consensus_${sample_id}_R1.fastq > consensus_${sample_id}_R1.fastq.gz;
    bgzip -c consensus_${sample_id}_R2.fastq > consensus_${sample_id}_R2.fastq.gz;
    """
}
process consensus_fastqQC{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/10_Consensus_fastqc_results", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    // maxForks 5

    input:
        tuple sample_id, "consensus_${sample_id}_R1.fastq.gz", "consensus_${sample_id}_R2.fastq.gz" from fastq_QC2 
    output:
        tuple sample_id, '*.{zip,html}' into fastqc_results
    script:
    """
    fastqc --quiet --outdir . consensus_${sample_id}_R1.fastq.gz consensus_${sample_id}_R2.fastq.gz
    """
}

process align_consensusFASTQ{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/11_post_consensus_alignment", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    maxForks 5
    input:
        tuple sample_id, "consensus_${sample_id}_R1.fastq.gz", "consensus_${sample_id}_R2.fastq.gz" from fastq_postconsensus
    output:
        tuple sample_id, "consensus_aligned_${sample_id}.sorted.bam" into markduplicates_ch, coverage_ch  
    script:
        """
        bwa mem -t ${params.CPU_BWA} -M -R \
            '@RG\\tID:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA\\tPM:MINISEQ\\tSM:${sample_id}' \\
            ${params.BWA_REPO}\
            consensus_${sample_id}_R1.fastq.gz consensus_${sample_id}_R2.fastq.gz > consensus_aligned_${sample_id}.sam;
        samtools view -S -b consensus_aligned_${sample_id}.sam | samtools sort -o consensus_aligned_${sample_id}.sorted.bam;
        
        """   
}

process MarkDuplicates{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/12_MarkDup", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    input:
        tuple sample_id, "consensus_aligned_${sample_id}.sorted.bam" from markduplicates_ch  
    output:
        tuple sample_id, "${sample_id}.dedup.bam", "${sample_id}.metrics.txt" into picard_CollectHsMetrics_ch, vardict, insert_metrics
    script:
    """
    samtools index consensus_aligned_${sample_id}.sorted.bam;
    picard MarkDuplicates -Xmx2G \
    INPUT=consensus_aligned_${sample_id}.sorted.bam \
    OUTPUT=${sample_id}.dedup.bam METRICS_FILE=${sample_id}.metrics.txt;
    """
}

process bedtools_coverage{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/13_bedtoolsCoverage", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    maxForks 1

    input:
        tuple sample_id, "consensus_aligned_${sample_id}.sorted.bam" from coverage_ch  
        file bedfile
    output:
        file "${sample_id}.cov_per_gene.txt" into bed_coverage
    script:
    """
    samtools index consensus_aligned_${sample_id}.sorted.bam;
    bedtools coverage -b consensus_aligned_${sample_id}.sorted.bam -a ${bedfile} > ${sample_id}.cov_per_gene.txt
    """
}

process CollectHsMetrics{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/14_HsMetrics", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    
    input:
        tuple sample_id, "${sample_id}.dedup.bam", "${sample_id}.metrics.txt" from picard_CollectHsMetrics_ch
        file bait_intervals
        file target_intervals    
    output:
        tuple sample_id, "${sample_id}.hs_metrics.txt", "${sample_id}.QC.txt" into hs_metrics
    script:
    """
    samtools index ${sample_id}.dedup.bam;
    picard CollectHsMetrics -Xmx2G \
    I=${sample_id}.dedup.bam \
    O=${sample_id}.hs_metrics.txt \
    R=${params.BWA_REPO}.fa \
    BAIT_INTERVALS=${bait_intervals} \
    TARGET_INTERVALS=${target_intervals}; 
    grep -B 1 "^${params.GENE_PANEL}" ${sample_id}.hs_metrics.txt > ${sample_id}.QC.txt
    """
}

process vardict{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/15_vardict_output", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    // maxForks 5 
    input:
        tuple sample_id, "${sample_id}.dedup.bam", "${sample_id}.metrics.txt" from vardict
        file bedfile
    output:
        tuple sample_id, "vardict_single_${sample_id}.vcf" into vardict_output
    script:
        """
        samtools index ${sample_id}.dedup.bam;
        vardict-java -t -G ${hg38} -q 20 -f ${params.AF_THR} -N $sample_id -b ${sample_id}.dedup.bam \
        -c 1 -S 2 -E 3 -g 4 ${bedfile} \
        | teststrandbias.R | var2vcf_valid.pl -N ${sample_id} -E -f ${params.AF_THR} > vardict_single_${sample_id}.vcf
        """
}

process annotate_vep{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/16_annotated_vcf", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    // maxForks 5
    input:
        tuple sample_id, "vardict_single_${sample_id}.vcf" from vardict_output
        file dir_cache 
        file fasta 
        file clinvar 
        file COSMIC
        file clinvar_index
        file COSMIC_INDEX
    output:
        tuple sample_id, "vardict_single_${sample_id}_annotated.vcf" into final_output
    script:
    """
    vep --stats_file ${sample_id}.vep_summary.html --cache\
    --dir_cache $dir_cache --refseq --force_overwrite --format vcf \
    --fasta $fasta --everything --pick   \
    --custom $clinvar,ClinVar,vcf,exact,0,ALLELEID,CLNSIG,GENEINFO,CLNREVSTAT,CLNDN,CLNHGVS \
    --custom $COSMIC,COSMIC,vcf,exact,0,GENE,STRAND,LEGACY_ID,SNP,CDS,HGVSC,HGVSG,CNT \
    -i vardict_single_${sample_id}.vcf --check_existing --symbol --vcf \
    -o vardict_single_${sample_id}_annotated.vcf --offline;
    """
}

process picard_CollectInsertSizeMetrics {
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/17_insert_metrics", mode: 'copy'
    // errorStrategy 'retry'
    // maxRetries 1
    // maxForks 5
    input:
        tuple sample_id, "${sample_id}.dedup.bam", "${sample_id}.metrics.txt" from insert_metrics

    output:
        tuple sample_id, "${sample_id}.insert_metrics.txt", "${sample_id}.histogram.pdf"

    script:
        """
        picard CollectInsertSizeMetrics \
        INPUT=${sample_id}.dedup.bam \
        OUTPUT=${sample_id}.insert_metrics.txt \
        HISTOGRAM_FILE=${sample_id}.histogram.pdf
        """
}


