#!/bin/bash -e
echo "Begin download TabAnno resource files under current directory..."
echo 

echo "Download reference file and its index:"
wget http://www.sph.umich.edu/csg/zhanxw/software/anno/resources/hs37d5.fa
wget http://www.sph.umich.edu/csg/zhanxw/software/anno/resources/hs37d5.fa.fai

echo "Download gene definition:"
wget http://www.sph.umich.edu/csg/zhanxw/software/anno/resources/refFlat_hg19.txt.gz

echo "Download TabAnno codon definition and annotation priority files:"
wget http://www.sph.umich.edu/csg/zhanxw/software/anno/codon.txt
wget http://www.sph.umich.edu/csg/zhanxw/software/anno/priority.txt

echo "Download finished."
echo 
echo "Example commands:"
PWD=`pwd`
echo "1. annotate VCF file and output to uncompressed VCF file"
echo "  [anno_executable] -i [input.vcf] -o [output.vcf] -r ${PWD}/hs37d5.fa -g ${PWD}/refFlat_hg19.txt.gz -p ${PWD}/priority.txt -c ${PWD}/codon.txt"
echo 
echo "2. annotate VCF file and output to compressed and indexed VCF file"
echo "  [anno_executable] -i [input.vcf] -o [output.vcf.gz] --indexOutput -r ${PWD}/hs37d5.fa -g ${PWD}/refFlat_hg19.txt.gz -p ${PWD}/priority.txt -c ${PWD}/codon.txt"
echo 
echo "3. annotate PLINK assoc files and output to plain text file"
echo "  [anno_executable] -i [input.assoc.assoc] -o [output.txt] --inputFormat plink -r ${PWD}/hs37d5.fa -g ${PWD}/refFlat_hg19.txt.gz -p ${PWD}/priority.txt -c ${PWD}/codon.txt"
