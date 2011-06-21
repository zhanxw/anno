/**
 * Main feature list
 * 1. support various annotion
 *    UPSTREAM / DOWNSTREAM 
 *    5PRIME_UTR / 3PRIME_UTR
 *    INTRON
 *    SYNONYMOUS / NONSYNONYMOUS
 *    STOP_GAIN / STOP_LOST
 *    START_GAIN / START_LOST ??
 *    SPLICE_SITE
 *    INTROGENETIC
 *    CODON_DELETION / CODON_INSERTION / FRAME_SHIFT
 *    
 * 2. input format are flexible
 *    support VCF, as well as user specified.
 * 3. output format are configurable.
 *    user specify
 * 4. output basis statisitcs:
 *    freq of each annotation  .anno.freq 
 *    freq of each base change .base.freq
 *    freq of codon change     .codon.freq
 *    freq of indel length     .indel.freq
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <string>
#include <vector>
#include <map>

#include "Argument.h"
#include "IO.h"
#include "StringUtil.h"

#include "Chromosome.h"

typedef enum AnnotationType{
    UPSTREAM = 0,
    DOWNSTREAM,
    UTR5,
    UTR3,
    INTRON,
    EXTRON,
    SYNONYMOUS,
    NONSYNONYMOUS,
    STOP_GAIN,
    STOP_LOST,
    START_GAIN,
    START_LOSE,
    NORMAL_SPLICE_SITE,
    ESSENTIAL_SPLICE_SITE
} AnnotationType;

char AnnotationString[][32]= {
    "Upstream", 
    "Downstream",
    "Utr5",
    "Utr3",
    "Intron",
    "Extron",
    "Synonymous",
    "Nonsynonymous",
    "Stop_Gain",
    "Stop_Lost",
    "Start_Gain",
    "Start_Lose",
    "Normal_Splice_Site",
    "Essential_Splice_Site"
};
// all are 1-based index, inclusive on boundaries.
struct Range{
    int start;
    int end;
    Range(int s, int e): start(s), end(e) {};
    Range(): start(-1), end(-1) {};
};
class Gene{
public:
    void readLine(const char* line) {
        std::vector< std::string > field;
        std::vector< std::string > exon_beg;
        std::vector< std::string > exon_end;
        int nf = stringTokenize(line, "\t", &field);
        if (nf < 11) {
            return;
        }
        this->name = field[0];
        this->chr = field[2];
        this->forwardStrand = (field[3] == "+" ? true: false);
        this->tx.start = toInt(field[4]);
        this->tx.end = toInt(field[5]);
        this->cds.start = toInt(field[6]);
        this->cds.end = toInt(field[7]);
        unsigned int nExon = toInt(field[8]);
        stringTokenize(field[9], ',', &exon_beg);
        stringTokenize(field[10], ',', &exon_end);
        for (unsigned int i = 0; i < nExon; i++ ){
            this->exon.push_back(Range(toInt(exon_beg[i]), toInt(exon_end[i])));
        }
    };
public:
    std::string name;
    std::string chr;
    bool forwardStrand;
    Range tx;
    Range cds;
    std::vector<Range> exon;
};
struct GeneAnnotationParam{
    GeneAnnotationParam():
        upstreamRange(500),
        downstreamRange(500),
        spliceIntoExon(2),
        spliceIntoIntron(2) {};
    int upstreamRange;      // upstream def
    int downstreamRange;    // downstream def
    int spliceIntoExon;     // essential splice site def
    int spliceIntoIntron;   // essentail splice site def
};

class GeneAnnotation{
public:
    void readGeneFile(const char* geneFileName){
        std::string line;
        std::vector<std::string> fields;
        LineReader lr(geneFileName);
        while (lr.readLine(&line) > 0) {
            stringStrip(&line);
            if (line.size()>0 and line[0] == '#' || line.size() == 0) continue; // skip headers and empty lines
            Gene g;
            g.readLine(line.c_str());
            this->geneList[g.chr].push_back(g);
        }
        return;
    }; 
    void openReferenceGenome(const char* referenceGenomeFileName) {

        // check if -bs.umfa file exists
        std::string umfaFileName = referenceGenomeFileName;
        stringSlice(&umfaFileName, 0, -3);
        umfaFileName += "-bs.umfa";
        FILE* fp = fopen(umfaFileName.c_str(), "r");
        if (fp == NULL) { // does not exist binary format, so try to create one
            fprintf(stdout, "Create binary reference genome file for the first run\n");
            this->gs.create();
        } else{
            fclose(fp);
        }
        this->gs.setReferenceName(referenceGenomeFileName);
        if (this->gs.open()) {
            fprintf(stderr, "Cannot open reference genome file %s\n", referenceGenomeFileName);
            exit(1);
        }
        if (!convert2Chromosome(&this->gs, &this->reference)) {
            fprintf(stderr, "Cannot use GenomeSequence by Chromosome.\n");
            exit(1);
        }
    };
    // we take a VCF input file for now
    void annotate(const char* inputFileName, const char* outputFileName){
        LineReader lr(inputFileName);
        std::vector<std::string> field;
        while (lr.readLineBySep(&field, "\t") > 0) {
            if (field.size() < 4) continue; 
        }
        return;
    };
    void setAnnotationParameter(GeneAnnotationParam& param) {
        this->param = param;
    };
private:
    // store results in start and end, index are 0-based, inclusive
    // find gene whose range plus downstream/upstream overlaps chr:pos 
    void findInRangeGene(unsigned int* start, unsigned int* end, const std::string& chr, unsigned int* pos) {
        return;
    };
    
    void annotateByGene(unsigned int geneIdx, const std::string& chr, const unsigned int& variantPos, const char* ref, const char* alt){
        Gene& g = this->geneList[chr][geneIdx];
        if (g.forwardStrand) {
            if (g.tx.start - param.upstreamRange < variantPos && variantPos < g.tx.start){
                this->annotation += AnnotationString[UPSTREAM];
            } else if (g.tx.end < variantPos && variantPos < g.tx.end + param.downstreamRange) {
                this->annotation += AnnotationString[DOWNSTREAM];
            } else {
            }
        } else { // backward strand
            
        }
        return;
    };
    /**@return -1: unknow type
     *          0: single mutation
     *          1: deletion
     *          2: insertion
     */
    int getMutationType(const char* ref, const char* alt) {
        unsigned int refLen = strlen(ref);
        unsigned int altLen = strlen(alt);
        if (refLen == altLen) 
            return 0;
        else if (refLen > altLen) 
            return 1;
        else if (refLen < altLen)
            return 2;
        return -1;
    };
    GeneAnnotationParam param;
    std::map <std::string, std::vector<Gene> > geneList;
    std::string annotation;
    GenomeSequence gs;
    std::map<std::string, Chromosome> reference;
};

int main(int argc, char *argv[])
{
    BEGIN_PARAMETER_LIST(pl)
        ADD_STRING_PARAMETER(pl, inputFile, "-i", "Specify input VCF file")
    END_PARAMETER_LIST(pl)
        ;
    
    pl.Read(argc, argv);
    pl.Status();
    
    GeneAnnotation ga;
    ga.readGeneFile("test_gene.txt");
    ga.openReferenceGenome("test.fa");
    ga.annotate("test.vcf", "test.output.vcf");

    return 0;
}
