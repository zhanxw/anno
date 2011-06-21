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
};
class Gene{
public:
    void ReadLine(const char* line);
public:
    std::string chr;
    Range tx;
    Range cds;
    std::vector<Range> exon;
    bool forwardStrand;
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
        return;
    }; 
    void openReferenceGenome(const char* referenceGenomeFileName) {
        this->gs.setReferenceName(referenceGenomeFileName);
        // check if -bs.umfa file exists
        std::string umfaFileName = referenceGenomeFileName;
        umfaFileName += "-bs.umfa";
        FILE* fp = fopen(umfaFileName.c_str(), "r");
        if (fp == NULL) { // does not exist binary format, so try to create one
            fprintf(stdout, "Create binary reference genome file for the first run\n");
            this->gs.create();
        } else{
            fclose(fp);
        }
        if (!this->gs.open()) {
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
    ga.readGeneFile("a");
    ga.openReferenceGenome("b");
    ga.annotate("input", "output");

    return 0;
}
