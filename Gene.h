#ifndef _GENE_H_
#define _GENE_H_
#include "Range.h"
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
        int cdsStart = toInt(field[6]) + 1;
        int cdsEnd = toInt(field[7]);
        if (this->length(cdsStart, cdsEnd) == 0) {
            this->isNonCodingGene = true;
        } else {
            this->isNonCodingGene = false;
        }
        unsigned int nExon = toInt(field[8]);
        stringTokenize(field[9], ',', &exon_beg);
        stringTokenize(field[10], ',', &exon_end);
        for (unsigned int i = 0; i < nExon; i++ ){
            this->exon.push_back(Range(toInt(exon_beg[i]) + 1, toInt(exon_end[i])));
        }
        if (!this->isNonCodingGene) {
            // we will assume this is forward strand
            // if this is not the case, we will swap utr5 and utr3

            if (name == "TMEM88B") {
                printf("\n");
            }
            
            // load till left of cdsBegin
            unsigned i = 0;
            for (i = 0; i < nExon; i++ ){
                int beg = exon[i].start;
                int end = exon[i].end;
                if (this->isInRange(cdsStart, beg, end)) {
                    if (beg != cdsStart) // avoid add empty range
                        this->utr5.push_back(Range(beg, cdsStart - 1));
                    break;
                } else {
                    this->utr5.push_back(Range(beg, end));
                }
            }
            // load cds region
            for (; i < nExon; i++) {
                int beg = exon[i].start;
                int end = exon[i].end;
                if (this->isInRange(cdsStart, beg, end)) {
                    beg = cdsStart;
                } 
                if (this->isInRange(cdsEnd, beg, end)) {
                    end = cdsEnd;
                    this->cds.push_back(Range(beg, cdsEnd));
                    break;
                } else {
                    this->cds.push_back(Range(beg,end));
                }
            }
            // load from cdsEnd to the end of exon
            for (; i < nExon; i++) {
                int beg = exon[i].start;
                int end = exon[i].end;
                if (this->isInRange(cdsEnd, beg, end)) {
                    if (cdsEnd != end) // avoid add empty range
                        this->utr3.push_back(Range(cdsEnd + 1, end));    
                } else {
                    this->utr3.push_back(Range(beg, end));    
                }
            }
            if (i != nExon) {
                assert (i == nExon);
            }
            if (!this->forwardStrand) {
                std::swap(this->utr5, this->utr3);
            }
        }

#if 0 
        // debug code 
        if (name == "DDX53") {
            assert ( 0 == getCDSLength() % 3 );
        }
#endif
#if 0
        // just for my curiosity
        if (!isNonCoding()) {
            if (name == "DDX53") {
                printf("\n");
            }
            printf("%s (%d) with %d exon has 5'utr (%d), 3'utr(%d), cds(%d), cds module 3(%d)\n",
                   name.c_str(),
                   getGeneLength(),
                   (int)exon.size(),
                   get5PrimeUTRLength(),
                   get3PrimeUTRLength(), getCDSLength(),
                   getCDSLength() % 3 );
        }
#endif
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
     * //@param utrPos will store the relative position of @param variousPos to the leftmost position of 5' UTR 
     * //@param utrLen will store the length of the 5'-UTR region
     */
    bool is5PrimeUtr(const int variantPos, int* utrPos, int* utrLen) {
        if (this->isNonCoding()) return false;
        if (this->isInRange(variantPos, this->utr5)){
            return true;
        };
        return false;
    };
    bool is3PrimeUtr(const int variantPos, int* utrPos, int* utrLen) {
        if (this->isNonCoding()) return false;
        if (this->isInRange(variantPos, this->utr3)){
            return true;
        };
        return false;
    };
    bool isExon(const int variantPos, int* exonNum){
        if (isNonCoding()) {
            if (this->isInRange(variantPos, this->exon))
                return true;
        } else {
            if (this->isInRange(variantPos, this->utr5) || 
                this->isInRange(variantPos, this->cds) || 
                this->isInRange(variantPos, this->utr3))
                return true;
        }
        return false;
    }
    int nextCodonPos(const int currentPos, int* cdsIdx, const int offset) {
        assert (offset == 1 || offset == -1);
        int nextPos = -1;
        if (offset == 1) {
            nextPos = currentPos + 1;
            if (!this->isInRange(nextPos, this->cds[*cdsIdx])) {
                (*cdsIdx) ++;
                if (*cdsIdx >= this->cds.size()) {
                    return -1;
                }
                nextPos = this->cds[*cdsIdx].start;
            } 
        } else {
            nextPos = currentPos - 1;
            if (!this->isInRange(nextPos, this->cds[*cdsIdx])) {
                (*cdsIdx) --;
                if (*cdsIdx < 0) {
                    return -1;
                }
                nextPos = this->cds[*cdsIdx].start;
            } 
        }            
        return nextPos;
    };
    /**
     * @return true: if codonPos[3] are all valid position
     * @param codonNum : which base (inclusive, 1-based) has mutation
     */
    bool calculatCodonPosition(const int variantPos, int* codonNum, int codonPos[3]){
        *codonNum = 0;
        if (this->forwardStrand) {
            unsigned int i;
            for (i = 0; i < this->cds.size() ; i++) {
                if (this->isInRange(variantPos, this->cds[i])){ 
                    *codonNum += variantPos - this->cds[i].start + 1;
                    break;
                } else {
                    *codonNum += this->cds[i].length();
                }
            }
            int n = (*codonNum) % 3;
            int cdsIdx = i;
            switch(n){
            case 0:
                codonPos[2] = variantPos;
                codonPos[1] = nextCodonPos(codonPos[2], &cdsIdx, -1);
                codonPos[0] = nextCodonPos(codonPos[1], &cdsIdx, -1);
                break;
            case 1:
                codonPos[0] = variantPos;
                codonPos[1] = nextCodonPos(codonPos[0], &cdsIdx, 1);
                codonPos[2] = nextCodonPos(codonPos[1], &cdsIdx, 1);
                break;
            case 2:
                codonPos[1] = variantPos;
                codonPos[2] = nextCodonPos(codonPos[1], &cdsIdx, 1);
                cdsIdx = i;
                codonPos[0] = nextCodonPos(codonPos[1], &cdsIdx, -1);
                break;
            }
        } else { // backward
            int i;
            for (i = this->cds.size() - 1; i >= 0 ; i--) {
                if (this->isInRange(variantPos, this->cds[i])){ 
                    *codonNum += this->cds[i].end - variantPos + 1;
                    break;
                } else {
                    *codonNum += this->cds[i].length();
                }
            }
            int n = (*codonNum) % 3;
            int cdsIdx = i;
            switch(n){
            case 0:
                codonPos[2] = variantPos;
                codonPos[1] = nextCodonPos(codonPos[2], &cdsIdx, +1);
                codonPos[0] = nextCodonPos(codonPos[1], &cdsIdx, +1);
                break;
            case 1:
                codonPos[0] = variantPos;
                codonPos[1] = nextCodonPos(codonPos[0], &cdsIdx, -1);
                codonPos[2] = nextCodonPos(codonPos[1], &cdsIdx, -1);
                break;
            case 2:
                codonPos[1] = variantPos;
                codonPos[2] = nextCodonPos(codonPos[1], &cdsIdx, -1);
                cdsIdx = i;
                codonPos[0] = nextCodonPos(codonPos[1], &cdsIdx, +1);
                break;
            }
        }
        if (codonPos[0] != -1 && codonPos[1] != -1 && codonPos[2] != -1)
            return true;
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
    int getTotalLength(const std::vector<Range>& v) {
        int l = 0;
        for (unsigned int i = 0; i < v.size() ; i++ )
            l += v[i].length();
        return l;
    };
    int getExonLength() {
        return this->getTotalLength(this->exon);
    };
    int getCDSLength() {
        return this->getTotalLength(this->cds);
    };
    int get5PrimeUTRLength() {
        return this->getTotalLength(this->utr5);
    };
    int get3PrimeUTRLength() {
        return this->getTotalLength(this->utr3);
    };
    int getGeneLength() {
        return this->tx.length();
    };
    bool isNonCoding() {
        return this->isNonCodingGene;
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
    bool isInRange(const int pos, const int beg, const int end) {
        if (beg > end) {
            fprintf(stdout, "in isInRange beg(%d) > end(%d).\n", beg, end);
        }
        if (beg <= pos && pos <= end) 
            return true;
        return false;
    };
    bool isInRange(const int pos, const Range& r) {
        return (this->isInRange(pos, r.start, r.end));
    };
    bool isInRange(const int pos, const std::vector<Range>& r) {
        for (unsigned int i = 0; i < r.size(); i++) {
            if (this->isInRange(pos, r[i]))
                return true;
        }
        return false;
    };
    /**
     * @return the total length from @param beg to @param end, inclusive on the boundaries
     */
    int length(int beg, int end) {
        if (beg > end+1) {
            fprintf(stdout, "in length beg(%d) > end(%d) + 1.\n", beg, end);
        }
        return (end - beg + 1);
    };
  public:
    std::string name;
    std::string chr;
    bool forwardStrand;
    Range tx;
    // used for nonCoding gene
    std::vector<Range> exon;
    // used for coding gene
    std::vector<Range> cds;
    std::vector<Range> utr5;
    std::vector<Range> utr3;
    bool isNonCodingGene;
};

#endif /* _GENE_H_ */