# rmats-nf: 

A pipeline to retrieve RNA sequencing data from SRA, perform alignment with hisat2 and characterize splicing events with rMATS 3.2.5.

# rmats-nf: Documentation

An overview of how the pipeline works, how to run it and a description of all of the different command-line flags.

## Usage

The typical command for running the pipeline is as follows:

```
nextflow run lifebit-ai/rmats-nf --accessionList 'accession_list.txt' -profile docker
```


### Mandatory arguments:
#### `--accessionList`

Path to input file with accession list to fetch from SRA

####  `--gencodeFile`

Path to input gencode GTF file

####  `--keyFile`

Path to a keyfile used to fetch restricted access datasets with SRAtools

#### `-profile`

Configuration profile to use. Can use multiple (comma separated)
Available: `short-test`, `key-test`

### Other arguments:

#### `--skipTrimming`

If used, trimming will be ommited.

#### `--readType`

Specifies what type of input reads are to be processed.
Accepted values: `single end`, `paired end`

#### `--readLength`

Specifies the read length
