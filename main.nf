/*
 * Copyright (c) 2019, Jackson Labs and the authors.
 *
 *   This file is part of 'rmats-nf' a pipeline repository to reproduce ... .
 *
 * Main rmats-NF pipeline script
 *
 * @authors
 * Anne Deslattes Mays
 * Pablo Prieto Barja <pablo.prieto.barja@gmail.com>
 */
log.info "RMATS - N F  ~  version 1.0"
log.info "====================================="
log.info "Accession list        : ${params.accessionList}"
log.info "Key file              : ${params.keyFile ? params.keyFile : 'Not provided'}"
log.info "Gencode annotation    : ${params.gencodeFile}"
log.info "Read type             : ${params.readType}"
log.info "Read length           : ${params.readLength}"
log.info "output                : ${params.output}"
log.info "\n"
def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run lifebit-ai/rmats-nf --accessionList 'accession_list.txt' -profile docker
    Mandatory arguments:
      --accessionList               Path to input file with accession list to fetch from SRA
      --gencodeFile                 Path to input gencode GTF file
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: short-test, key-test, ...
    Generic:
      --readType                    Specifies what type of input reads are to be processed: single end or paired end
      --readLength                  Specifies the read length
    """.stripIndent()
}

/*********************************
 *      CHANNELS SETUP           *
 *********************************/
key_file = file(params.keyFile)
Channel
    .fromPath( params.accessionList )
    .splitText()
    .map{ it.trim() }
    .dump(tag: 'AccessionList content')
    .set { accessionIDs }
Channel
        .value(file(params.gencodeFile))
        .ifEmpty { error "Cannot find any Gencode GTF annotaton for parameter --gencode: ${params.gencodeFile}" }
        .set { gencodeFile }
if (!params.readLength) {
        exit 1, "Cannot find read length value assigned"
}
if (!params.readType) {
        exit 1, "Cannot find read type value assigned"
}

**********************************
 *      PIPELINE PROCESSES        *
 **********************************/
/*
 * Get accession samples
 */
process getAccession {
    tag "${accession}"
    input:
    val accession from accessionIDs
    file keyFile from key_file
    output:
    set val(accession), file("*.fastq.gz") into readFiles
    script:
    def vdbConfigCmd = keyFile.name != 'NO_FILE' ? "vdb-config --import ${keyFile} ./" : ''
    """
    $vdbConfigCmd
    fasterq-dump $accession --threads ${task.cpus} --split-3
    pigz *.fastq
    """
}


process trimming {
    tag "${accession}"
    input:
    set val(accession), file(reads) from readFiles
    output:
    set val(accession), file("*_trimmed.fastq.gz") into trimmedFiles
    script:
    if (params.readType == "single") {
      """
      fastp -w ${task.cpus} -i ${reads} -b ${params.readLength} -o ${accession}_trimmed.fastq.gz
      """
    } else {
      """
      fastp -w ${task.cpus} -i ${reads[0]} -I ${reads[1]} -b ${params.readLength} -B  ${params.readLength} -o ${accession}_R1_trimmed.fastq.gz -O ${accession}_R2_trimmed.fastq.gz
      """
    }
}

/*
 * Mapping
 */
process mapping {
    publishDir = [path: "${params.output}/alignments", mode: 'copy', overwrite: 'true' ]
    tag "mapping: $reads"
    cache 'lenient'
    input:
    set val(name), file(reads) from trimmedFiles
    output:
    set val(name), file("${name}.bam") into hisat2Bams
    file("${name}.hisat2_summary.txt") into hisat2Multiqc
    script:
    if (params.readType == "single") {
                """
            hisat2 -x ${params.hisat2index} \
                   -U ${reads} \
                   -p ${task.cpus} \
                   --met-stderr \
                   --new-summary \
                   --summary-file ${name}.hisat2_summary.txt \
                   | samtools view -bS - > ${name}.bam
                """
    } else {
        """
            hisat2 -x ${params.hisat2index} \
                   -1 ${reads[0]} \
                   -2 ${reads[1]} \
                   -p ${task.cpus} \
                   --met-stderr \
                   --new-summary \
                   --summary-file ${name}.hisat2_summary.txt \
                   | samtools view -bS - > ${name}.bam
        """
    }
}

/*
 * SortBam
 */
process sortbam {
    publishDir = [path: "${params.output}/sorted_alignments", mode: 'copy', overwrite: 'true' ]
    tag "sortbam: $name"
    input:
    set val(name), file(bam) from hisat2Bams
    output:
    set val(name), file("${name}.sorted.bam") into sortedBams
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

*
 * MarkDuplicates
 */
process markduplicates {
    publishDir = [path: "${params.output}/marked_duplicates", mode: 'copy', overwrite: 'true' ]
    tag "markdups: $name"
    input:
    set val(name), file(bam) from sortedBams
    output:
    set val(name), file("${name}.sorted.nodup.bam") into markedDups, bamstatsSamples
    file("${name}.metric.txt") into markduplicatesMultiqc
    script:
    // Runtime java Mark duplicates options
    markdup_java_options="-Xmx30g"
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
    """
}

markedDups
        .flatten()
        .collate( 4 )
        .set { pairedSamples }
/*
 * Generate BAM statistics
 */
process bamstats {
    publishDir = [path: "${params.output}/flagstats", mode: 'copy', overwrite: 'true' ]
    tag "bamstats: $name"
    input:
    set val(name), file(bam) from bamstatsSamples
    output:
    file("${name}.sorted.nodup.bam.flagstat.txt") into flagstatMultiqc
    script:
    """
    samtools flagstat $bam > ${name}.sorted.nodup.bam.flagstat.txt
    samtools view $bam | cut -f 10 | perl -ne 'chomp;print length(\$_) . "\n"' | sort | uniq -c >> ${name}.sorted.nodup.bam.flagstat.txt
    """
}

process paired_rmats {
    publishDir "${params.output}/paired_rmats", mode: 'copy',
        saveAs: {filename ->
              if (filename.indexOf("fromGTF.novelEvents") > 0) null
              else if (filename.indexOf("fromGTF") > 0) "${sample1Name}_${sample2Name}/${filename}"
              else null
        }
    tag "paired_rmats: ${sample1Name}_${sample2Name}"
    input:
    set val(sample1Name), file(sample1Bam), val(sample2Name), file(sample2Bam) from pairedSamples
    file(gencodeGtf) from gencodeFile
    output:
    file 'fromGTF.*.txt' into splicingEvents
    script:
    """
    ls $sample1Bam > b1.txt
    ls $sample2Bam > b2.txt
    rmats.py --nthread ${task.cpus} --b1 b1.txt --b2 b2.txt --gtf $gencodeGtf --od ./ -t ${params.readType} --readLength ${params.readLength} --statoff
    """
}
/*
    sampleCountsSave.sh ./ ${sample1Name} ${sample2Name}
*/
