#!/bin/bash
fqgz=$1
fstem=${fqgz%_R1_001.fastq.gz}
bcf=$2
zcat $fqgz | fastx_clipper -a AGATCGGAAGAGCACACGTCTGAACTC -ncv \
	| fastx_barcode_splitter.pl --eol --bcfile $bcf --prefix ${fstem}_ --suffix "_disamb.fastq"
