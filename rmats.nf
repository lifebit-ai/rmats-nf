/*
 * Copyright (c) 2019, Jackson Labs and the authors.
 *
 *   This file is part of 'rmats-nf'.
 *
 * Main rmats-NF pipeline script
 *
 * @authors
 * Anne Deslattes Mays
 * Pablo Prieto Barja <pablo.prieto.barja@gmail.com>
 */


params.genome_file = "$baseDir/example/genome/chrX_reduced.fa"
params.input_folder  = "$baseDir/example/reads"
params.reads_extension = "*_{1,2}.fastq.gz"
params.sample_id     = "ERR188383"
params.output        = "results/"


log.info "RMATS - N F  ~  version 1.0"
log.info "====================================="
log.info "Sample ID             : ${params.sample_id}"
log.info "Input folder          : ${params.input_folder}/${params.reads_extension}"
log.info "Genome                : ${params.genome_file}"
log.info "output                : ${params.output}"
log.info "\n"

/*********************************
 *      CHANNELS SETUP           *
 *********************************/

Channel
    .fromPath ( params.genome_file )
    .set { genomes }



Channel
    .fromFilePairs( "${params.input_folder}/${params.reads_extension}", size: -1)
    .ifEmpty { error "Cannot find any reads matching: ${params.input_folder}/${params.reads_extension}"}
    .set { read_files } 



/**********************************
 *      PIPELINE PROCESSES        *
 **********************************/


/*
 * Indexing
 */

process genome_index {
    publishDir = [path: "${params.output}/index", mode: 'copy', overwrite: 'true' ]

    input:
    file genome_file from genomes

    output:
    set val ("genome_index"), file("genomeindex") into genome_index

    script:
    """
    hisat2-build ${genome_file} genome_index
    mkdir genomeindex
    mv genome_index* genomeindex/.
    """
}

/*
 * Mapping
 */

process mapping {
    publishDir = [path: "${params.output}/mapped_sams", mode: 'copy', overwrite: 'true' ]
    tag "mapping: $reads"

    input:
    set val(index_name), file(index_dir) from genome_index.first()
    set val(name), file(reads) from read_files

    output:
    set val(name), file("${name}.bam") into hisat2_bams
    file("${name}.hisat2_summary.txt") into hisat2_multiqc

    script:
    """
    hisat2 -x $index_dir/genome_index \\
           -1 ${reads[0]} \\
           -2 ${reads[1]} \\
           -p ${task.cpus} \\
           --met-stderr \\
           --new-summary \\
           --summary-file ${name}.hisat2_summary.txt \\
           | samtools view -bS - > ${name}.bam
    """
}

/*
 * SortBam
 */

process sortbam {
    publishDir = [path: "${params.output}/mapped_bams", mode: 'copy', overwrite: 'true' ]
    tag "sortbam: $name"

    input:
    set val(name), file(bam) from hisat2_bams

    output:
    set val(name), file("${name}.sorted.bam") into sorted_bams

    script:
    avail_mem=""
    """   
    samtools sort \\
        $bam \\
        -@ ${task.cpus} ${avail_mem} \\
        -o ${name}.sorted.bam

    samtools index ${name}.sorted.bam
    """
}

/*
 * MarkDuplicates
 */

process markduplicates {
    publishDir = [path: "${params.output}/mapped_bams", mode: 'copy', overwrite: 'true' ]
    tag "markdups: $name"

    input:
    set val(name), file(bam) from sorted_bams

    output:
    set val(name), file("${name}.sorted.nodup.bam") into marked_dups
    file("${name}.metric.txt") into markduplicates_multiqc

    script:
    // Runtime java Mark duplicates options
    markdup_java_options="-Xms2g -Xmx4g"
    """   

    picard ${markdup_java_options} MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${name}.sorted.nodup.bam \\
        METRICS_FILE=${name}.metric.txt \\
        REMOVE_DUPLICATES=true \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT
    
    samtools index ${name}.sorted.nodup.bam
    """
}

/*
 * Final touch
 */

process bamstats {
    publishDir = [path: "${params.output}/flagstats", mode: 'copy', overwrite: 'true' ]
    tag "bamstats: $name"

    input:
    set val(name), file(bam) from marked_dups

    output:
    file("${name}.sorted.nodup.bam.flagstat.txt") into flagstat_multiqc

    script:
    """   
    samtools flagstat $bam > ${name}.sorted.nodup.bam.flagstat.txt
    samtools view $bam | cut -f 10 | perl -ne 'chomp;print length(\$_) . "\n"' | sort | uniq -c >> ${name}.sorted.nodup.bam.flagstat.txt
    """
}

/*
 * MultiQC
 */

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
       
    input:
    file('*') from hisat2_multiqc.collect()
    file('*') from markduplicates_multiqc.collect()
    file('*') from flagstat_multiqc.collect()
    
    output:
    file("*") into viz
     
    script:
    """
    multiqc . 
    """
}

workflow.onComplete {
        println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.output/multiqc_report.html\n" : "Oops .. something went wrong" )
}
