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
log.info "Genome                : ${params.genomeFile}"
log.info "HISAT2 index          : ${params.hisat2_index}"
log.info "Read type             : ${params.readType}"
log.info "Read length        	: ${params.readLength}"
log.info "output                : ${params.output}"
log.info "\n"

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run lifebit-ai/rmats-nf --accessionList --genomeFile 'genome.fa'  -profile docker
    Mandatory arguments:
      --accessionList               Path to input file with accession list to fetch from SRA
      --genomeFile                  Path to input genome FASTA file to index and map against
      --gencodeFile                 Path to input gencode GTF file

      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: short-test, key-test, ...
    Generic:
      --readType                    Specifies what type of input reads are to be processed: single end or paired end
      --readLength                  Specifies the read length
    References:                     If not specified in the configuration file or you wish to overwrite any of the references.
      --genome                      Name of iGenomes reference
      --hisat2_index                Path to HiSAT2 index : avoids re-building it
    """.stripIndent()
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*********************************
 *      CHANNELS SETUP           *
 *********************************/

keyFile = file(params.keyFile)

Channel
    .fromPath( params.accessionList )
    .splitText()
    .map{ it.trim() }
    .dump(tag: 'AccessionList content')
    .set { accessionIDs }

Channel
    .fromPath ( params.genomeFile )
    .set { genome }

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

if (params.hisat2_index) {
  if (hasExtension(params.hisat2_index, 'gz')) {
    hs2_indices_gz = Channel
        .fromPath("${params.hisat2_index}", checkIfExists: true)
        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
  } else {
    hs2_indices = Channel
        .fromPath("${params.hisat2_index}*", checkIfExists: true)
        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
  }
}

/**********************************
 *      PIPELINE PROCESSES        *
 **********************************/

if (params.hisat2_index ) {
	process gunzip_hisat_index {
	    tag "$gz"

	    input:
	    file gz from hs2_indices_gz

	    output:
	    file "**/*.ht2*" into hs2_indices

	    script:
	    // Use tar as the hisat2 indices are a folder, not a file
	    """
	    tar -xzvf ${gz}
	    """
	}
}



/*
 * Indexing
 */
if (!params.hisat2_index ) {
	process genome_index {
	    publishDir = [path: "${params.output}/index", mode: 'copy', overwrite: 'true' ]

	    input:
	    file genomeFile from genome

	    output:
	    file "${genomeFile.baseName}.*.ht2*" into hs2_indices

	    script:
	    """
	    hisat2-build -p ${task.cpus} ${genomeFile} ${genomeFile.baseName}.hisat2_index
	    """
	}
}


/*
 * Get accession samples
 */

process getAccession {
    
    tag "${accession}"

    input:
    val accession from accessionIDs

    output:
    set val(accession), file("*.fastq.gz") into readFiles

    script:
    def vdbConfigCmd = keyFile.name != 'NO_FILE' ? "vdb-config --import ${keyFile} ./" : ''
    """
    $vdbConfigCmd
    #prefetch -X ${params.spotsNumber} -N ${task.cpus} $accession
    #fastq-dump -N ${task.cpus} -I --origfmt --split-3 $accession --gzip
    fasterq-dump $accession --threads ${task.cpus} --split-3
    gzip *.fastq
    """
}

/*
 * Mapping
 */

process mapping {
    publishDir = [path: "${params.output}/alignments", mode: 'copy', overwrite: 'true' ]
    tag "mapping: $reads"

    input:
    file(hs2_indices) from hs2_indices.collect()
    set val(name), file(reads) from readFiles

    output:
    set val(name), file("${name}.bam") into hisat2Bams
    file("${name}.hisat2_summary.txt") into hisat2Multiqc

    script:
    index_base = hs2_indices[0].toString() - ~/.\d.ht2l?/

    if (params.readType == "single") {
		"""
	    hisat2 -x $index_base \
	           -U ${reads} \
	           -p ${task.cpus} \
	           --met-stderr \
	           --new-summary \
	           --summary-file ${name}.hisat2_summary.txt \
	           | samtools view -bS - > ${name}.bam
		"""
    } else {
    	"""
	    hisat2 -x $index_base \
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

/*
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
    markdup_java_options="-Xms3g -Xmx7g"
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


/*
 *
 */

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
    file("rmats_output/fromGTF.*.txt") into splicingEvents

    script:
    """   
    rmats.py --nthread ${task.cpus} --b1 ${sample1Bam} --b2 ${sample2Bam} --gtf $gencodeGtf --od rmats_output -t ${params.readType} --readLength ${params.readLength} --statoff
    sampleCountsSave.sh rmats_output ${sample1Name} ${sample2Name}
    """
}



/*
 * 
    Add missing commands to create all the matrices in the previous command maybe

process create_matrices {

    publishDir "${params.output}/paired_rmats", mode: 'copy', 
        saveAs: {filename ->
              if (filename.indexOf("fromGTF.novelEvents") > 0) null
              else if (filename.indexOf("fromGTF") > 0) "${sample1Name}_${sample2Name}/${filename}"
              else null
        }
    tag "paired_rmats: ${sample1Name}_${sample2Name}"

    input:

    output:

    script:
    """
    create_matrices_from_files.sh $astype $jctype $cnttype $inputDir $splitNum
    rm rmats_matrix*
    """
}
 */

 // Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}
