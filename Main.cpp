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
#include <string>
#include <vector>
#include <map>

class FileReader{
public:
    FileReader(const char* fileName):
        fp(NULL) {
        this->fp = fopen(fileName, "r");
        if (!this->fp) {
            fprintf(stderr, "ERROR: Cannot open %s\n", fileName);
        }
    };
    ~FileReader() {
        if (this->fp) {
            fclose(fp);
        }
    };
    // return number of characters read.
    unsigned int readLine(std::string* line) {
        
    };
    // return number of fields read.
    unsigned int readLineBySep(std::vector<std::string>* fields) {
    };
    // check eof 
    bool isEof() {
        return (feof(this->fp) != 0);
    }

private:
    FILE* fp;
};
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

char AnnotationString[][16]= {
    "UPSTREAM", 
    "DOWNSTREAM"
};
// all are 1-based index, inclusive on boundaries.
struct Range{
    int start;
    int end;
};
class Gene{
public:
    void ReadLine(const char* line);
private:
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
    void readGeneFile(const char* geneFileName){
        return;
    }; 
    void annotate(const char* inputFileName, const char* outputFileName){
        return;
    };
private:
    // store results in start and end, index are 0-based, inclusive
    // find gene whose range plus downstream/upstream overlaps chr:pos 
    void findInRangeGene(unsigned int* start, unsigned int* end, const char* chr, unsigned int* pos) {
        return;
    };
    void annotateByGene(unsigned int geneIdx, const unsigned int& variantPos, const char* ref, const char* alt){
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
};

int main(int argc, char *argv[])
{
    return 0;
}
