#!/bin/bash
#
# 2020 January 20
#
# Purpose:  rMATS.4.0.2 now outputs variable length files for the results from the aligning
#           of RNAseq files to the junctions.
#
#           This causes problems with merging the results together into a single very useful matrix.
#
# Resolution: Normalize these to the orgin - definition file.
#           rMATS scans the transcriptome definition file, novel or reference (aka gencode.v32.annotation.gtf)
#           and creates 10 definition files:
#
#           fromGTF.A3SS.txt
#           fromGTF.A5SS.txt
#           fromGTF.MXE.txt
#           fromGTF.RI.txt
#           fromGTF.SE.txt
#           fromGTF.novelEvents.A3SS.txt
#           fromGTF.novelEvents.A5SS.txt
#           fromGTF.novelEvents.MXE.txt
#           fromGTF.novelEvents.RI.txt
#           fromGTF.novelEvents.SE.txt
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
echo "astype = $astype"
echo "jctype = $jctype"
echo "cnttype= $cnttype"

#
# definitions
#
all="*"
joiner="."
normalized="n"
ending="txt"
filelist="filelist"
fromGTF="fromGTF"
id="ID"
# convert the provided lower case as type to upper
upperastype=$( echo "$astype" | tr '[a-z]' '[A-Z]')
filelist=$all$joiner$astype$joiner$jctype$joiner$cnttype$joiner$ending


fromGTF_file=$fromGTF$joiner$upperastype$joiner$ending
fromGTF_id=$fromGTF$joiner$upperastype$joiner$id$joiner$ending

echo "fromGTF file  = $fromGTF_file"
echo "fromGTF ID file  = $fromGTF_id"

#
# make the ID file
#
cut -f 1 $fromGTF_file > $fromGTF_id

echo "filelist     = $filelist"
echo "fromGTF_id   = $fromGTF_id"

files=$(ls $filelist)

for file in $files;do
    normalized_file=$(echo $file| sed s'/.txt/.n.txt/g');
#    echo "file is $file";
#    echo "normalized file is $normalized_file";
    awk 'NR==FNR{A[$1]=$1;B[$1]=$2;next} {if (A[$1]!="") result=B[$1];else result=0; print $1" "result}' $file $fromGTF_id > $normalized_file;
done


