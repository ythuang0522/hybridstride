# Introduction
The HybridStriDe Assembler aims to correct, overlap, and assembly low-quality PacBio reads in conjunction with high-quality Illumina reads. The entire implementation is derived from Illumina-only [StriDe][3] assembler, which was based on Simpson's [SGA][1] and Li's [ropebwt2][2]. The HybridStriDe is still frequently updated. This repository serves as a quick download version for the submission. 

# Executable version
A precompiled version under Linux 64bit (StriDe_Linux_64bit) can be directly downloaded and executed. 

# Compile by yourself
To compile StriDe assembler in your specific environment, type 

      1. ./autogen.sh 
      2. ./configure
      3. make

An executable program called stride will be found under the StriDe folder.

# Script Execution
A sample script is provided for all-in-one execution. Please modify the mandatory paramters (PacBio Coverage, Illumina covergae and read length) before execution. Given one PacBio reads (PB.fa) and one Illumina paired-end reads (R1.fq, R2.fq), type

	bash HybridStriDe.sh R1.fq R2.fq PB.fa

# Step by step exceution
If the all-in-one script does not fit your datasets, you can execute each command line in the sample bash script and adjust the parameters for your own needs.

HybridStride preprocesss --discard-quality -p 1 R1.fq R2.fq -o reads.fa
HybridStride index -t reads.fa
HybridStride index -t PB.fa
...


[1]: https://github.com/jts/sga
[2]: https://github.com/lh3/ropebwt2
[3]: https://github.com/ythuang0522/StriDe
