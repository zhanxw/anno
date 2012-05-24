= Introduction = Anno is short for "annotation", it is used for annotate
variants. Our goal is to provide abundant information for the variant.
For example, it can annotate different transcripts of the same genes.
Anno support various file format: VCF file, plain file, plink
association output file and epacts file.

= Where to Find It =

The compiled executable file is at:
/net/fantasia/home/zhanxw/anno/executable/anno

The source code is located at:/net/fantasia/home/zhanxw/anno You can
type 'make release' to compile your own executable file. Type "make
test1" or "make test2" will demonstrate the command line to annotate
example VCF files.

= Usage =

== Command line == After you obtain the anno executable (either by
compiling the source code or by downloading the pre-compiled binary
file), you will find the executable file under executable/anno.

Here is the anno help page by invoking anno without any command line
arguments:

some\_linux\_host \> executable/anno
..............................................\
 ... G(ene) A(nnotation) ...\
 ... Xiaowei Zhan, Goncalo Abecasis ...\
 ... zhanxw@umich.edu ...\
 ... Jun 2011 ... ................................................

Required Parameters -i : Specify input VCF file -o : Specify output VCF
file -g : Specify gene file Optional Parameters -r : Specify reference
genome position -f : Specify gene file format (default: refFlat, other
optio ns knownGene) -c : Specify codon file (default: codon.txt) -u :
Specify upstream range (default: 50) -d : Specify downstream range
(default: 50) --se : Specify splice into exon range (default: 3) --si :
Specify splice into intron range (default: 8) Please specify input file

== Input files ==

anno runs on the input VCF file specified on the command-line using flag
'-i'.

Additionally, you need to specify gene file using flag '-g'. You can use
the default refFlat file (using HG19 genome build):
/net/fantasia/home/zhanxw/anno/refFlat\_hg19.txt.gz

== Parameters ==

Some of the command line parameters are described here, but most are
self explanatory.

\*Reference genome file

-r Specify a FASTA format reference genome file.

Specify ''-r'' option enable anno to give more detailed information, for
example, instead of annotating a variant as exon, it will tell you which
codon, which exon that variant locates, whether its
synonymous/non-synonymous and etc.

Anno requires Fasta index file and that will save running memory and
speed up annotation. You can use "samtools faidx XXX.fa" to generate
Fasta index file.

For example, you can specify Fasta file of the whole genome and use "-r
/data/local/ref/karma.ref/human.g1k.v37.fa"

\*Gene file format

Currently, anno support gene file in refFlat format. A prepared list of
all gene obtained of UCSC website is:
"/net/fantasia/home/zhanxw/anno/refFlat\_hg19.txt.gz" . To use that
file, you use flag ''-g''.

As anno support refFlat file by default, you can use refFlat format
without specify gene format flag ''-f''.

To use knownGene or refGene format, you need to specify both ''-g'' and
''-f'' flat to tell anno which gene file and which format it is. For
example, ''-g /net/fantasia/home/zhanxw/anno/knownGene.txt.gz -f
knownGene''.

\*Codon file

Codon file can tell the relationship between triplet codons and amino
acids. A default file is located in:
''/net/fantasia/home/zhanxw/anno/codon.txt''. If you have special codon
file, you can specify that using flag ''-c'', otherwise, anno will use
the default codon file:

''default codon file'' \# DNA codon table. Stop codon is coded as: Stp,
O, Stop \#Codon AminoAcid Letter FullName AAA Lys K Lysine AAC Asn N
Asparagine AAG Lys K Lysine AAT Asn N Asparagine ACA Thr T Threonine ...
\*Annotation ranges -u how far apart from 5'-end of the gene are counted
as upstream -d how far apart from 3'-end of the gene are counted as
upstream -se how far apart from 5'-end of the gene are counted as
upstream -si how far apart from 3'-end of the gene are counted as
upstream

= Example =

anno can annotate a VCF file and also output statistics of 4 frequency
table: annotation type; base change; codon change; indel size. More
details will be given below.

== Built-in example ==

In example/ folder, you can see test.vcf, which is a toy example. You
can invoke anno using the following command line:

cd example; ./anno -i test.vcf -r test.fa -g test.gene.txt -c
../codon.txt -o test.out.vcf

Sample outputs are listed below:

1)  Annotated VCF file, ''test.out.vcf''

\#VCF\_test \#from
http://www.sph.umich.edu/csg/liyanmin/vcfCodingSnps/Tutorial.shtml
\#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA12891 NA12892 NA12878 1
3 . A G 50 depth=20
ANNO=GENE1/CODING\_GENE:+:Exon:Utr5:Normal\_Splice\_Site|GENE3/CODING\_GENE:-:Exon:Utr3:Normal\_Splice\_Site|GENE2/NON\_CODING\_GENE:+:Upstream
GT:GQ:GD 1/0:31:12 0/0:28:14 1 5 . C A 50 depth=20
ANNO=GENE1/CODING\_GENE:+:Exon:Nonsynonymous(CCT/Pro/P-\>CAT/His/H:Base3/30:Codon1/10:Exon1/5):Normal\_Splice\_Site|GENE3/CODING\_GENE:-:Exon:Nonsynonymous(AGG/Arg/R-\>ATG/Met/M:Base30/30:Codon10/10:Exon5/5):Normal\_Splice\_Site|GENE2/NON\_CODING\_GENE:+:Upstream
GT:GQ:GD 1/0:31:12 0/0:28:14 ...

The annotation results are stored in the INFO column after ANN= tag. The
annotation format is defined as following:

"|" separates different transcripts, e.g. in the first line, chromosome
1 position 3, there are 3 annotations:
"GENE1/CODING\_GENE:+:Exon:Utr5:Normal\_Splice\_Site" and
"GENE3/CODING\_GENE:-:Exon:Utr3:Normal\_Splice\_Site" and
"GENE2/NON\_CODING\_GENE:+:Upstream"

":" separates within gene annotation in the following order: gene,
strand, exon/intron, details.

2)  Statistics files:

Four frequency table will be generated after annotation. For example:

''test.out.vcf.anno.frq''\
 Stop\_Loss 1 Utr5 2 Utr3 2 CodonRegion 2 CodonGain 2 Frameshift 2
Synonymous 3 StructuralVariation 3 Noncoding 3 Nonsynonymous 4 Deletion
6 Upstream 6 Insertion 6 Essential\_Splice\_Site 8 Downstream 8 Intron
12 Exon 21 Normal\_Splice\_Site 25

''test.out.vcf.base.frq'' A-\>G 1 T-\>C 1 T-\>G 2 A-\>C 2 C-\>A 5

''test.out.vcf.codon.frq'' Arg-\>Met 1 Pro-\>Thr 1 Arg-\>Arg 1 Pro-\>His
1 Gly-\>Gly 1 Pro-\>Pro 1 Stp-\>Tyr 1 Leu-\>Val 1

''test.out.vcf.indel.frq'' 1 1 -4 1 3 1 -3 1

= Contact =

Questions and requests should be sent to Xiaowei Zhan
([mailto:zhanxw@umich.edu zhanxw@umich.edu]) or Goncalo Abecasis
([mailto:goncalo@umich.edu goncalo@umich.edu])
