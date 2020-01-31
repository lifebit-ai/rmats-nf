#!/bin/bash
#
# 2018 March 20
# 2020 January 21
#
# Update -- need to use .n.txt - the normalized files to guarantee all input files are the same length
#           due to the fact that rMATS.4.0.2 outputs variable data output formats.
#
# Purpose:  To collate all the single run samples into a single matrix.
#           The call includes the specification of which one to build
#           se, a3ss, a5ss, mxe or ri
#           Input will be a bucket containing the single column sample matrix by the junction IDs
#           for each of the different splicing events
#
# Assumption: The bucket contains the rMATS RNASeq-MATS.py output
#
# Inputs:
#
# rMATS RNASeq-MATS.py produces 10 different output types which get assembled into as type junction ID by sample ID matrices
#
# Alternative Splice Site Types are: (se, a3ss, a5ss, mxe, ri)
# This is input as ARGV1 into variable 'astype'
#  a.  Skipped Exon events (se),
#  b.  Alternative 3' splice site (a3ss),
#  c.  Alternative 5' splice site (a5ss),
#  d.  Mutually exclusive exon (mxe),
#  e.  and retention intron (ri)
#
#
# There are two different kinds of junction counts
# This is input as ARGV2 into variable 'jctype'
#  a. jc = junction counts - reads that cross the junction
#  b. jcec = junction counts plus reads on the target (such as included exon
#
# And the count type -- there are 5 types
# This is input as ARGV into variable 'cnttype'
#  a. inclusion levels (percent spliced in)
#  b. included junction counts (ijc)
#  c. skipped junction counts (sjc)
#  d. inclusion length (inclen)
#  e. skipped length (skiplen)
#
#------------------------------------------------------------------------------------------------------
astype=$1
jctype=$2
cnttype=$3
#inputDir=$4
splitNum=$4

echo "astype   = $astype"
echo "jctype   = $jctype"
echo "cnttype  = $cnttype"
#echo "inputDir = $inputDir"
echo "splitNum = $splitNum"
#
# definitions
#
joiner="."
ending="n.txt"
astypeSpecificFile=$astype$joiner$jctype$joiner$cnttype$joiner$ending
inputFiles=*$joiner$astypeSpecificFile
splitFiles=$astypeSpecificFile$joiner
rmats_matrix="rmats_matrix"
rmats_final="rmats_final"
filelist="filelist"
create_matrix_script="create-matrix.py"
create_matrix=$(which "$create_matrix_script" 2>/dev/null)
[ $? -gt 0 -a -f "$0" ] && create_matrix="$create_matrix_script"

rmats_filelist=$rmats_final$joiner$astypeSpecificFile$joiner$filelist$joiner$ending
rmats_final_matrix=$rmats_final$joiner$astypeSpecificFile

echo "inputFiles         = $inputFiles"
echo "astypeSpecificFile = $astypeSpecificFile"
echo "splitFiles         = $splitFiles"
echo "rmats_matrix       = $rmats_matrix"

ls $inputFiles > $astypeSpecificFile

split -l $splitNum $astypeSpecificFile $splitFiles

files=$(ls $splitFiles*)
echo "files = $files"
echo "create-matrix is $create_matrix"

for file in $files; do rmats_name="$rmats_matrix.$file"; echo $rmats_name; python $create_matrix $file $rmats_name ' '; done

rmats_matrices=$(ls rmats_matrix*)
echo "rmats_matrices     = $rmats_matrices"
echo "rmats_filelist     = $rmats_filelist"
echo "rmats_final_matrix = $rmats_final_matrix"
echo $rmats_matrices | tr ' ' '\n' > $rmats_filelist

python $create_matrix $rmats_filelist $rmats_final_matrix ','

