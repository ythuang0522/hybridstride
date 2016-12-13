#!/bin/bash

#### Mandatory parameter parts #####
# Path of binary stride executable file
SD=./HybridStride_Linux64bit

# number of threads
threads=24

# coverage of PacBio (C) and Illumina (c), e.g., 30x and 100x respectively.
PBCoverage=30
IlluCoveragec=100

# Illumina read length
r=100

#### The main script parts ####
# Illumina 1st and 2nd ends of paired reads
R1=$1
R2=$2

# PacBio reads
PB=$3

# minimum overlap parameter
ovl=749

startTime=`date "+%s"`

# preprocess short reads into interleaved format:
$SD preprocess --discard-quality -p 1 $1 $2 -o reads.fa

# Construct FM-index of Illumina reads
$SD index -a ropebwt2 -t $threads reads.fa

# Optinally, the user may correct the Illumina reads first or not, depending on the time you want to spend. The following script assumes error correction of illumina is performed.
$SD correct -a overlap -t $threads -k 31 -x 3 reads.fa -o READ.ECOLr.fasta

# Construct FM-index from corrected Illumina reads
$SD index -t $threads READ.ECOLr.fasta

# Hybrid error correction for PacBio reads using Illumina reads
$SD pbhc -p READ.ECOLr -f PB -t $threads -c $IlluCoverage -C $PBCoverage -r $r PB.fa
mv PB.PBHybridCor.fa merged.fa

# Filter redundant reads
$SD index -a ropebwt2 -t $threads merged.fa
$SD filter -t $threads merged.fa

# Overlap Computation stage
$SD overlap -m $ovl -e 0.012 -l 50 merged.filter.pass.fa -t $threads

# Assembly stage
# -i is the median or N50 PacBio read length
$SD asmlong merged.filter.pass.asqg.gz

endTime=`date "+%s"`
elapsed=$((endTime-startTime))
printf "Total time : %i:%02i:%02i\n" $((elapsed/3600)) $(((elapsed/60)%60)) $(($elapsed%60))
