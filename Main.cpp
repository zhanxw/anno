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

/**
 *           @
 * 0-base: 0 1 2 3 4 5 
 * 1-base: 1 2 3 4 5 6
 *             ^
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
    EXON,
    SYNONYMOUS,
    NONSYNONYMOUS,
    STOP_GAIN,
    STOP_LOST,
    START_GAIN,
    START_LOSE,
    NORMAL_SPLICE_SITE,
    ESSENTIAL_SPLICE_SITE,
    INTROGENIC
} AnnotationType;

char AnnotationString[][32]= {
    "Upstream", 
    "Downstream",
    "Utr5",
    "Utr3",
    "Intron",
    "Exon",
    "Synonymous",
    "Nonsynonymous",
    "Stop_Gain",
    "Stop_Lost",
    "Start_Gain",
    "Start_Lose",
    "Normal_Splice_Site",
    "Essential_Splice_Site",
    "Introgenic"
};
// all are 1-based index, inclusive on boundaries.
// NOTE: UCSC use 0-based start and 1-based end, this is convenient to calculate length
//       but for complex case, this causes confusion.
struct Range{
    int start;
    int end;
    Range(int s, int e): start(s), end(e) {};
    Range(): start(-1), end(-1) {};
    int getLength() { return (end - start + 1); };
};
class Gene{
public:
    void readLine(const char* line) {
        std::vector< std::string > field;
        std::vector< std::string > exon_beg;
        std::vector< std::string > exon_end;
        int nf = stringTokenize(line, "\t", &field);
        if (nf < 11) { 
            static int nTimeError = 0;
            fprintf(stderr, "Unable to read this gene from: %s\n", line);
            if (nTimeError++ > 10) {
                fprintf(stderr, "Too many errors, now quiting...\n");
                exit(1);
            }
            return;
        }
        this->name = field[0];
        this->chr = field[2];
        this->forwardStrand = (field[3] == "+" ? true: false);
        this->tx.start = toInt(field[4]) + 1;
        this->tx.end = toInt(field[5]);
        this->cds.start = toInt(field[6]) + 1;
        this->cds.end = toInt(field[7]);
        unsigned int nExon = toInt(field[8]);
        stringTokenize(field[9], ',', &exon_beg);
        stringTokenize(field[10], ',', &exon_end);
        for (unsigned int i = 0; i < nExon; i++ ){
            this->exon.push_back(Range(toInt(exon_beg[i]) + 1, toInt(exon_end[i])));
        }
    };
    /**
     *@return true if @param pos is in upstream and return how far it is from the beginning of the gene
     */
    bool isUpstream(const int pos, const int upstreamRange, int* dist) {
        if (this->forwardStrand) {
            if (this->tx.start - upstreamRange < pos && pos < this->tx.start) {
                *dist = this->tx.start - pos;
                return true;
            }
        } else {
            if (this->tx.end < pos && pos < this->tx.end + upstreamRange) {
                *dist = this->tx.end - pos;
                return true;
            }
        }
        return false;
    };
    bool isDownstream(const int pos, const int downstreamRange, int* dist) {
        if (this->forwardStrand) {
            if (this->tx.end < pos && pos < this->tx.end + downstreamRange) {
                *dist = pos - this->tx.end;
                return true;
            }
        } else {
            if (this->tx.start - downstreamRange < pos && pos < this->tx.start) {
                *dist = this->tx.start - pos;
                return true;
            }
        }
        return false;
    };
    /**
     * @return true is @param variousPos is i 5'-UTR region, 
     * @param utrPos will store the relative position of @param variousPos to the leftmost position of 5' UTR 
     * @param utrLen will store the length of the 5'-UTR region
     */
    bool is5PrimeUtr(const int variantPos, int* utrPos, int* utrLen) {
        if (this->isNonCoding()) return false;
        if (this->forwardStrand) {
            if (this->exon[0].start <= variantPos && variantPos < this->cds.start) {
                *utrPos = variantPos - this->exon[0].start;
                *utrLen = this->cds.start - this->exon[0].start + 1;
                return true;
            }
        } else {// backward strand
            if (this->cds.end < variantPos && variantPos <= this->exon[this->exon.size() - 1].end) {
                *utrPos = variantPos - (this->cds.end + 1); // +1: since this->cds.end is not in 5'-UTR region
                *utrLen = this->exon[this->exon.size() - 1].end - this->cds.end + 1;
                return true;
            }
        }
        return false;
    };
    bool is3PrimeUtr(const int variantPos, int* utrPos, int* utrLen) {
        if (this->isNonCoding()) return false;
        if (this->forwardStrand) {
            if (this->cds.end < variantPos && variantPos <= this->exon[this->exon.size() - 1].end) {
                *utrPos = variantPos - (this->cds.end + 1) ; // +1: since this->cds.end is not in 5'-UTR region
                *utrLen = this->exon[this->exon.size() - 1].end - this->cds.end + 1;
                return true;
            }
        } else {// backward strand
            if (this->exon[0].start <= variantPos && variantPos < this->cds.start) {
                *utrPos = variantPos - this->cds.start - variantPos;
                *utrLen = this->cds.start - this->exon[0].start + 1;
                return true;
            }
        }
        return false;
    };
    bool isExon(const int variantPos, int* exonNum){
        if (isNonCoding()) return false;
        *exonNum = 0;
        for (unsigned int i = 0; i < this->exon.size() ; i++) {
            if (this->isInRange(variantPos, this->exon[i].start, this->exon[i].end)) {
                *exonNum = i;
                return true;
            }
        }
        return false;
    }
    bool calculatCodon(const int variantPos, const int exonNum, int* codonNum, int codonPos[3]){
        // calculate which codon it hits
        // we did not count by 3, we count total bases from the start of CDS to variantPos
        *codonNum = 0;
        if (this->forwardStrand) {
            if (exonNum == 0) {
                *codonNum += this->length(this->cds.start, variantPos);
            } 
            for (unsigned int i = 1; i < exonNum; i++) {
                *codonNum += this->length(this->exon[i].start, this->exon[i].end);
            }
            *codonNum += this->length(this->exon[exonNum].start, variantPos);
        } else {
            if (exonNum == this->exon.size() - 1) {
                *codonNum += this->length(variantPos, this->exon[this->exon.size() - 1].end);
            } 
            for (unsigned int i = this->exon.size() - 2; i != exonNum; i--) {
                *codonNum += this->length(this->exon[i].start, this->exon[i].end);
            }
            *codonNum += this->length(variantPos, this->exon[exonNum].end);
        }
        return false;
    };
    bool isIntron(const int variantPos, int* intronNum){
        // strand is not an issue here
        for (unsigned int i = 1; i < this->exon.size(); i++) {
            if (this->exon[i-1].end < variantPos && this->exon[i].start) {
                return true;
            }
        }
        return false;
    };
    bool isSpliceSite(const int variantPos, int spliceIntoExon, int spliceIntoIntron, bool* isEssentialSpliceSite){
        *isEssentialSpliceSite = false;
        unsigned int exonNumber = this->exon.size();
        // first check splice into exon
        if (this->isInRange(variantPos, this->exon[0].end - (spliceIntoExon - 1), this->exon[0].end)) {
            return true;
        } 
        if (this->isInRange(variantPos, this->exon[exonNumber - 1].start, this->exon[exonNumber - 1].start + (spliceIntoExon - 1))){
                return true;
        }
        for (unsigned int i = 1; i < exonNumber - 1; i ++) {
            if (this->isInRange(variantPos, this->exon[i].start, this->exon[i].start + (spliceIntoExon - 1)))
                return true;
            if (this->isInRange(variantPos, this->exon[i].end - (spliceIntoExon - 1), this->exon[i].end))
                return true;
        }
        // check splice into intron (also mark isEssentialSpliceSite)
        // we define essential splice site is (GU...AG) in the intron, and next to exon
        // and GU, AG both have length 2.
        for (unsigned int i = 0; i < exonNumber - 1; i++ ) {
            if (this->isInRange(variantPos, this->exon[i].end+1, this->exon[i].end+1+2)) {
                *isEssentialSpliceSite = true;
                return true;
            } else if (this->isInRange(variantPos, this->exon[i+1].start - 1 - 2, this->exon[i+1].start - 1)) {
                *isEssentialSpliceSite = true;
                return true;
            }
            if (this->isInRange(variantPos, this->exon[i].end+1, this->exon[i].end + 1 + (spliceIntoIntron - 1))) {
                return true;
            } else if (this->isInRange(variantPos, this->exon[i+1].start - 1 - (spliceIntoIntron - 1), this->exon[i+1].start - 1)) {
                return true;
            }
        }
        return false;
    };
    int getTotalExonLength() {
        int l = 0;
        for (int i = 0; i < this->exon.size(); i++) {
            l += this->exon[i].getLength();
        }
    };
    bool isNonCoding() {
        if (this->cds.start == this->cds.end) {
            return true;
        }
        return false;
    };
    /**
     * given @param variantPos, @param exonNum (which exon the variantPos lies),
     * @param codonNum (how many bases from the begining of cds to 
     */
    void calculateCodonPos(int exonNum, int codonNum, int variantPos, int codonPos[3]) {
        
    };
    /**
     * @return true if @param pos is in the range [@param beg, @param end] (inclusive on the boundaries).
     */
    bool isInRange(int pos, int beg, int end) {
        if (beg > end) {
            fprintf(stdout, "in isInRange beg(%d) > end(%d).\n", beg, end);
        }
        if (beg <= pos && pos <= end) 
            return true;
        return false;
    };
    /**
     * @return the total length from @param beg to @param end, inclusive on the boundaries
     */
    int length(int beg, int end) {
        if (beg > end) {
            fprintf(stdout, "in length beg(%d) > end(%d).\n", beg, end);
        }
        return (end - beg + 1);
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
        spliceIntoExon(3),
        spliceIntoIntron(8) {};
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
        this->gs.setReferenceName(referenceGenomeFileName);
        if (fp == NULL) { // does not exist binary format, so try to create one
            fprintf(stdout, "Create binary reference genome file for the first run\n");
            this->gs.create();
        } else{
            fclose(fp);
        }
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
        std::string annotationString;
        LineReader lr(inputFileName);
        std::vector<std::string> field;
        while (lr.readLineBySep(&field, "\t") > 0) {
            if (field.size() < 4) continue;
            std::string chr = field[0];
            int pos = toInt(field[1]);
            std::string ref = field[3];
            std::string alt = field[4];
            int geneBegin;
            int geneEnd;
            this->findInRangeGene(&geneBegin, &geneEnd, field[0], &pos);
            if (geneEnd < 0) continue;
            annotationString.clear();
            for (unsigned int i = geneBegin; i <= geneEnd; i++) {
                this->annotateByGene(i, chr, pos, ref, alt, &annotationString);
                printf("%s %d ref=%s alt=%s has annotation: %s\n",
                       chr.c_str(), pos, 
                       ref.c_str(), alt.c_str(),
                       annotationString.c_str());
            }
        }
        return;
    };
    void setAnnotationParameter(GeneAnnotationParam& param) {
        this->param = param;
    };
private:
    // store results in start and end, index are 0-based, inclusive
    // find gene whose range plus downstream/upstream overlaps chr:pos 
    void findInRangeGene(int* start, int* end, const std::string& chr, int* pos) {
        *start = -1;
        *end = -1;
        std::vector<Gene>& g = this->geneList[chr];
        if (g.size() == 0) {
            return;
        } 
        unsigned int gLen = g.size();
        for (unsigned int i = 0; i < gLen ; i++) {
            if (g[i].forwardStrand) { // forward strand
                if (g[i].tx.start - param.upstreamRange < (*pos) && (*pos) < g[i].tx.end + param.downstreamRange) {
                    if (*start < 0) {
                        *start = i;
                    }
                    *end = i;
                }
            } else { // reverse strand
                if (g[i].tx.start - param.downstreamRange < *pos && *pos < g[i].tx.end + param.upstreamRange) {
                    if (*start < 0) {
                        *start = i;
                    }
                    *end = i;
                }
            }
        }
        // printf("start = %d, end = %d \n", *start, *end);
        return;
    };
    
    void annotateByGene(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt,
                        std::string* annotationString){
        Gene& g = this->geneList[chr][geneIdx];
        std::vector<std::string> annotation;
        int dist2Gene;
        int utrPos, utrLen;
        int exonNum; // which exon
        int codonNum; // which codon
        int codonPos[3] = {0, 0, 0}; // the codon position
        std::string codonRef;
        std::string codonAlt;
        AnnotationType type; // could be one of 
                             //    SYNONYMOUS / NONSYNONYMOUS
                             //    STOP_GAIN / STOP_LOST
                             //    START_GAIN / START_LOST ??
        int intronNum; // which intron
        bool isEssentialSpliceSite;
        if (g.isUpstream(variantPos, param.upstreamRange, &dist2Gene)) {
            annotation.push_back(AnnotationString[DOWNSTREAM]);
        } else if (g.isDownstream(variantPos, param.upstreamRange, &dist2Gene)) {
            annotation.push_back(AnnotationString[DOWNSTREAM]);
        } else if (g.isExon(variantPos, &exonNum)){//, &codonNum, codonPos)) {
            annotation.push_back(AnnotationString[EXON]);
            if (g.is5PrimeUtr(variantPos, &utrPos, &utrLen)) {
                annotation.push_back(AnnotationString[UTR5]);
            } else if (g.is3PrimeUtr(variantPos, &utrPos, &utrLen)) {
                annotation.push_back(AnnotationString[UTR3]);
            } else if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                if (isEssentialSpliceSite)
                    annotation.push_back(AnnotationString[ESSENTIAL_SPLICE_SITE]);
                else 
                    annotation.push_back(AnnotationString[NORMAL_SPLICE_SITE]);
            }
        } else if (g.isIntron(variantPos, &intronNum)) {
            annotation.push_back(AnnotationString[INTRON]);
            if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                if (isEssentialSpliceSite)
                    annotation.push_back(AnnotationString[ESSENTIAL_SPLICE_SITE]);
                else 
                    annotation.push_back(AnnotationString[NORMAL_SPLICE_SITE]);
            }
        } else {
            annotation.push_back("Intergenic");
        }
        annotationString->assign(stringJoin(annotation, ":"));
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
    ga.readGeneFile("test.gene.txt");
    ga.openReferenceGenome("test.fa");
    ga.annotate("test.vcf", "test.output.vcf");

    return 0;
}
