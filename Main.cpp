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

class Chromosome
class Sequence{
public:
    void open(const char* fileName) {
        
    }
private:
    std::map < std::string, std::string > sequence;
    
};
class FileReader{
public:
    FileReader(const char* fileName):
        fp(NULL) {
        this->open(fileName);
    };
    ~FileReader() {
        if (this->fp) {
            fclose(fp);
        }
    };
    
    // return number of characters read.
    // when reading an empty line, will return 1, as we read '\n', however, line will be empty
    // when reading the end, we will return 0
    unsigned int readLine(std::string* line) {
        if (this->isEof()) return 0;
        assert(line);
        line->clear();
        char c;
        unsigned nRead = 0;
        while (true) {
            c = this->getc();
            if (c == EOF) {
                return nRead;
            } else if (c == '\r') {
                // skip this
                continue;
            } else if (c == '\n') {
                ++nRead;
                return nRead;
            } else { // normal characters
                ++nRead;
                line->push_back(c);
            }
        }   
        assert(false); // should not reach here
        return 0;
    };
    // return number of fields read.
    // when reading an empty line, will return 1, meaning 1 field are read, although its content is empty
    // when reading to the EOF, will 0. 
    unsigned int readLineBySep(std::vector<std::string>* fields, const char* seq) {
        if (this->isEof()) return 0;
        assert(fields);
        assert(seq);
        fields->clear();
        char c;
        std::string s;
        while (true) {
            c = this->getc();
            if (c == EOF) {
                fields->push_back(s);
                return fields->size();
            } else if (c == '\r') {
                // skip this
                continue;
            } else if (c == '\n') {
                fields->push_back(s);
                return fields->size();
            } else if (strchr(seq, c) != NULL) { // separator
                fields->push_back(s);
                s.clear();
            } else { // normal characters
                s.push_back(c);
            }
        }   
        assert(false); // should not reach here
        return 0;
    };
    // get a char, if EOF, return EOF
    int getc(){
        return ::getc(this->fp);
    }
    // check eof 
    bool isEof() {
        return (feof(this->fp) != 0);
    }
    // open
    FILE* open(const char* fileName) {
        this->fp = fopen(fileName, "r");
        if (!this->fp) {
            fprintf(stderr, "ERROR: Cannot open %s\n", fileName);
        }
    }
    // close 
    void close() {
        if (this->fp) {
            fclose(fp);
        }
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
    void readGeneFile(const char* geneFileName){
        return;
    }; 
    void openReferenceGenome(const char* referenceGenomeFileName) {
        this->gs.setReferenceName(referenceGenomeFileName);
        if (!this->gs.open()) {
            fpritnf(stderr, "Cannot open reference genome file %s\n", referenceGenomeFileName);
            exit(1);
        }
    };
    // we take a VCF input file for now
    void annotate(const char* inputFileName, const char* outputFileName){
        FileReader fr(inputFileName);
        std::vector<std::string> field;
        while (fr.readLineBySep(&field, "\t") > 0) {
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
            } else if (g.tx.end < variantPos && variantPos < g.tx.end + downstreamRange) {
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
};

int main(int argc, char *argv[])
{
    // FileReader fr("Makefile");
    // // std::string line;
    // // while (fr.readLine(&line) > 0) {
    // //     fprintf(stdout, "%s\n", line.c_str());
    // // }
    // fr.close();
    
    // fr.open("Makefile");
    // std::vector<std::string> f;
    // while(fr.readLineBySep(&f, "\t") > 0) {
    //     for (unsigned int  i = 0; i < f.size(); i++) {
    //         if (i)
    //             printf("\t");
    //         printf("%s", f[i].c_str());
    //     }
    //     printf("\n");
    // }

    return 0;
}
