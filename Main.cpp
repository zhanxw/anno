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
 * Example, 1-index range 3-4, inclusive on the boundary
 * in UCSC setting, this range is coded as 2-4 (see below)
 *             @
 * 0-base: 0 1 2 3 4 5 
 * 1-base: 1 2 3 4 5 6
 *               ^
 * the tricky thing is 0-length range, in such case,
 * the UCSC coded the start and end as the same value
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "Argument.h"
#include "IO.h"
#include "StringUtil.h"

#include "Gene.h"
#include "GeneFormat.h"

char complementBase(const char c){
    switch (c){
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    case 'a':
        return 't';
    case 't':
        return 'a';
    case 'g':
        return 'c';
    case 'c':
        return 'g';
    default:
        return 'N';
    }
};

void reverseComplementTriplet(char s[3]){
    char t = s[0];
    s[0] = complementBase(s[2]);
    s[2] = complementBase(t);
    s[1] = complementBase(s[1]);
};

class Codon{
public:
    bool open(const char* codonFile) {
        LineReader lr(codonFile);
        std::string line;
        std::vector<std::string> f;
        while(lr.readLine(&line)>0) {
            if (line.size() > 0 && line[0] != '#') {
                stringTokenize(line, '\t', &f);
                this->codon2aa[f[0]] = f[1];
                this->codon2letter[f[0]] = f[2];
                this->codon2fullName[f[0]] = f[3];
            };
        }
    };
    const std::string& toAA(const char s[3]) {
        std::string key;
        key.push_back(s[0]);
        key.push_back(s[1]);
        key.push_back(s[2]);
        return safeAccess(this->codon2aa, key, Codon::unknownAA);
    };
public:
    static std::string unknownAA;
    static std::string unknownLetter;
    static std::string unknownFullName;
private:
    const std::string& safeAccess( const std::map<std::string, std::string>& data, 
                                   const std::string& key,
                                   const std::string& defaultValue) const {
        std::map<std::string, std::string>::const_iterator it = data.find(key);
        if (it == data.end())
            return defaultValue;
        return it->second;
    }

private:
    std::map<std::string, std::string> codon2aa;        // three letter amino acid
    std::map<std::string, std::string> codon2letter;    // amino acio letter
    std::map<std::string, std::string> codon2fullName;  // full amino acid name
    
};
std::string Codon::unknownAA = "N/A";
std::string Codon::unknownLetter = "*";
std::string Codon::unknownFullName ="UnknownAminoAcid";

class GenomeSequence{
public:
    /**
     * @return true: if loads successful
     */
    bool open(const char* fileName){
        LineReader lr(fileName);
        std::string line;
        std::string chr = "NoName";
        while(lr.readLine(&line) > 0) {
            if (line.size() > 0) {
                if (line[0] == '>') {
                    // new chromosome
                    unsigned int idx = line.find_first_of(' ');
                    chr = line.substr(1, idx - 1);
                    this->data[chr] = "";
                } else {
                    stringStrip(&line);
                    this->data[chr] += line;
                }
            }
        }
        return true;
    };
    int size() const {
        return this->data.size();
    };
    std::string& getChromosome(const char* c){
        return this->data[c];
    };
    const std::string& operator[] (const std::string& c) {
        return this->data[c];
    }
    bool exists(const std::string& c){
        if (this->data.find(c) != this->data.end())
         return true;
        return false;
    }
public:
    std::map<std::string, std::string> data;
};

typedef enum AnnotationType{
    UPSTREAM = 0,
    DOWNSTREAM,
    UTR5,
    UTR3,
    INTRON,
    EXON,
    SNV,                    /*SNV contain the following 6 types, it appears when there is no reference.*/
    SYNONYMOUS,
    NONSYNONYMOUS,
    STOP_GAIN,
    STOP_LOSS,
    START_GAIN,
    START_LOSS,
    NORMAL_SPLICE_SITE,
    ESSENTIAL_SPLICE_SITE,
    INTROGENIC
} AnnotationType;

const char* AnnotationString[]= {
    "Upstream", 
    "Downstream",
    "Utr5",
    "Utr3",
    "Intron",
    "Exon",
    "SNV",
    "Synonymous",
    "Nonsynonymous",
    "Stop_Gain",
    "Stop_Loss",
    "Start_Gain",
    "Start_Loss",
    "Normal_Splice_Site",
    "Essential_Splice_Site",
    "Introgenic"
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
    void setFormat(const std::string& f) {
        if (f == "refFlat") {
            this->format.setRefFlatFormat();
        } else if (f == "knownGene") {
            this->format.setUCSCKnownGeneFormat();
        } else {
            fprintf(stderr, "Unknown format!\nNow quitting...\n");
        }
    };

    void openGeneFile(const char* geneFileName){
        std::string line;
        std::vector<std::string> fields;
        LineReader lr(geneFileName);
        while (lr.readLine(&line) > 0) {
            stringStrip(&line);
            if (line.size()>0 && line[0] == '#' || line.size() == 0) continue; // skip headers and empty lines
            Gene g;
            g.readLine(line.c_str(), this->format);
            this->geneList[g.chr].push_back(g);
        }
        // make sure genes are ordered
        this->sortGene();
        return;
    }; 
    void openCodonFile(const char* codonFileName) {
        this->codon.open(codonFileName);
    };
    void openReferenceGenome(const char* referenceGenomeFileName) {
        this->gs.open(referenceGenomeFileName);
        return;
    };
    // we take a VCF input file for now
    void annotate(const char* inputFileName, const char* outputFileName){
        std::string annotationString;
        LineReader lr(inputFileName);
        std::vector<std::string> field;
        std::string line;
        while(lr.readLine(&line) > 0) {
            if (line.size() == 0 || line[0] == '#') continue;
            stringTokenize(line, "\t", &field);
            if (field.size() < 4) continue;
            std::string chr = field[0];
            int pos = toInt(field[1]);
            std::string ref = toUpper(field[3]);
            std::string alt = toUpper(field[4]);
            std::vector<unsigned> potentialGeneIdx;
            this->findInRangeGene(field[0], &pos, &potentialGeneIdx);
            if (potentialGeneIdx.size() == 0) continue;
            annotationString.clear();
            for (unsigned int i = 0; i < potentialGeneIdx.size(); i++) {
                this->annotateByGene(potentialGeneIdx[i], chr, pos, ref, alt, &annotationString);
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
    // make sure genes are ordered
    void sortGene() {
        std::map<std::string, std::vector<Gene> >:: iterator it;
        for (it = this->geneList.begin(); it != this->geneList.end(); it ++){
            std::sort( it->second.begin(), it->second.end(), GeneStartCompareLess);
        }
    };
    // store results in @param potentialGeneIdx
    // find gene whose range plus downstream/upstream overlaps chr:pos 
    void findInRangeGene(const std::string& chr, int* pos, std::vector<unsigned int>* potentialGeneIdx) {
        assert(potentialGeneIdx);
        potentialGeneIdx->clear();

        std::vector<Gene>& g = this->geneList[chr];
        unsigned int gLen = g.size();
        if (gLen == 0) {
            return;
        } 
        int maxDist = (param.upstreamRange > param.downstreamRange) ? param.upstreamRange : param.downstreamRange;
        Range r ((*pos - maxDist), (*pos + maxDist));
        for (unsigned int i = 0; i < gLen; i++ ){
            if (g[i].tx.start <= r.start) {
                if (g[i].tx.end < r.start){
                    continue;
                } else 
                    potentialGeneIdx->push_back(i);
            } else if (r.isInRange(g[i].tx.start)) {
                potentialGeneIdx->push_back(i);
            } else {
                break;
            }
        }
#if 0
        for (unsigned int i = 0 ; i < potentialGeneIdx->size() ; i++){
            printf("%d, ", (*potentialGeneIdx)[i]);
        }
        printf("\n");
#endif
        return;
    };
    // void annotateCodon(const std::string& chr, const int variantPos, const int codonNum, const int codonPos[3], 
    //                    const std::string& ref, const std::string& alt,
    //                    char refTriplet[3], std::string *refAAName,
    //                    char altTriplet[3], std::string *altAAName,
    //                    AnnotationType *at) {

    //     // obtain AA name
    //     *refAAName = this->codon.toAA(refTriplet);
    //     *altAAName = this->codon.toAA(altTriplet);

    //     *at = deter(*refAAName, *altAAName, codonNum);
    //     return ;
    // };
    /**
     * fill the actual base in @param refTriplet and @param altTriplet
     * we consider @param forwardStrand, so for forward strand, we copy from this->reference, 
     * or, we copy the reverse complement from this->reference
     */
    void fillTriplet(const std::string& chr, const int variantPos, const int codonPos[3], bool forwardStrand,  
                     const std::string& ref, const std::string& alt, 
                     char refTriplet[3], char altTriplet[3]) {
        assert(ref.size() == 1 && alt.size() == 1);
        const std::string& seq = this->gs[chr];
        if (codonPos[0] < 0 || codonPos[2] > seq.size()) {
            refTriplet[0] = refTriplet[1] = refTriplet[2] = 'N';
            altTriplet[0] = altTriplet[1] = altTriplet[2] = 'N';
        } else {
            refTriplet[0] = seq[codonPos[0] - 1];
            refTriplet[1] = seq[codonPos[1] - 1];
            refTriplet[2] = seq[codonPos[2] - 1];
            altTriplet[0] = (variantPos != codonPos[0]) ? seq[codonPos[0] - 1] : alt[0];
            altTriplet[1] = (variantPos != codonPos[1]) ? seq[codonPos[1] - 1] : alt[0];
            altTriplet[2] = (variantPos != codonPos[2]) ? seq[codonPos[2] - 1] : alt[0];
        }
        if (!forwardStrand) {
            reverseComplementTriplet(refTriplet);
            reverseComplementTriplet(altTriplet);
        }
    };
    AnnotationType determineSNVType(const std::string& refAAName, const std::string& altAAName, const int codonNum){
        if (refAAName == Codon::unknownAA || altAAName == Codon::unknownAA) {
            return SNV;
        } else if (refAAName == "Stp" && altAAName != "Stp") {
            return STOP_LOSS;
        } else if (refAAName != "Stp" && altAAName == "Stp") {
            return STOP_GAIN;
        } else if (refAAName == "Met" && altAAName != "Met" && codonNum <= 3) {
            return START_LOSS;
        } else if (refAAName != "Met" && altAAName == "Met" && codonNum <= 3) {
            return START_GAIN;
        } else if (refAAName == altAAName) {
            return SYNONYMOUS;
        } else {
            return NONSYNONYMOUS;
        }
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
            if (!g.isNonCoding()) {
                if (g.is5PrimeUtr(variantPos, &utrPos, &utrLen)) {
                    annotation.push_back(AnnotationString[UTR5]);
                } else if (g.is3PrimeUtr(variantPos, &utrPos, &utrLen)) {
                    annotation.push_back(AnnotationString[UTR3]);
                } else { // cds part has base change
                    if (g.calculatCodonPosition(variantPos, &codonNum, codonPos)) {
                        char refTriplet[3];
                        char altTriplet[3];
                        std::string refAAName;
                        std::string altAAName;
                        AnnotationType annotationType;
                        // when reference genome is provided
                        if (this->gs.size() >0){
                            // this->annotateCodon(chr, variantPos, codonNum, codonPos, ref, alt,
                            //                     refTriplet, &refAAName,
                            //                     altTriplet, &altAAName,
                            //                     &annotationType);
                            this->fillTriplet(chr, variantPos, codonPos, g.forwardStrand, ref, alt, refTriplet, altTriplet);

                            refAAName = this->codon.toAA(refTriplet);
                            altAAName = this->codon.toAA(altTriplet);
                            annotationType = this->determineSNVType(refAAName, altAAName, codonNum);

                            std::string s;
                            s = AnnotationString[annotationType];
                            s += refTriplet[0];
                            s += refTriplet[1];
                            s += refTriplet[2];
                            s += ":";
                            s += refAAName;
                            s += "->";
                            s += altTriplet[0];
                            s += altTriplet[1];
                            s += altTriplet[2];
                            s += ":";
                            s += altAAName;
                            annotation.push_back(s);
                        } else {
                            annotation.push_back(AnnotationString[SNV]);
                        }
                    } else {
                        annotation.push_back(AnnotationString[SNV]);
                    }
                }
            }
            // check splice site
            if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                if (isEssentialSpliceSite)
                    annotation.push_back(AnnotationString[ESSENTIAL_SPLICE_SITE]);
                else 
                    annotation.push_back(AnnotationString[NORMAL_SPLICE_SITE]);
            }
        } else if (g.isIntron(variantPos, &intronNum)) {
            annotation.push_back(AnnotationString[INTRON]);
            // check splice site
            if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                if (isEssentialSpliceSite)
                    annotation.push_back(AnnotationString[ESSENTIAL_SPLICE_SITE]);
                else 
                    annotation.push_back(AnnotationString[NORMAL_SPLICE_SITE]);
            }
        } else {
            //annotation.push_back("Intergenic");
        }
        *annotationString = g.name;
        annotationString->push_back(':');
        *annotationString+=(stringJoin(annotation, ":"));
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
    Codon codon;
    GeneFormat format;
};
int main(int argc, char *argv[])
{

    BEGIN_PARAMETER_LIST(pl)
        ADD_STRING_PARAMETER(pl, inputFile, "-i", "Specify input VCF file")
        ADD_STRING_PARAMETER(pl, geneFile, "-g", "Specify gene file")
        ADD_STRING_PARAMETER(pl, referenceFile, "-r", "Specify reference genome position")
        ADD_STRING_PARAMETER(pl, geneFileFormat, "-f", "Specify gene file format (default: refFlat, other options knownGene)")
    END_PARAMETER_LIST(pl)
        ;
    
    pl.Read(argc, argv);
    if (FLAG_geneFileFormat.size() == 0) {
        FLAG_geneFileFormat = "refFlat";
    }        
#if 1
    GeneAnnotation ga;
    if (FLAG_inputFile.size() == 0) {
        pl.Help();
        fprintf(stderr, "Please specify input file\n");
        exit(1);
    }
    if (FLAG_geneFile.size() == 0) {
        pl.Help();
        fprintf(stderr, "Please specify gene file\n");
        exit(1);
    }

    pl.Status();
    if (FLAG_referenceFile.size() != 0) {
        ga.openReferenceGenome(FLAG_referenceFile.c_str());
        ga.openCodonFile("codon.txt");
    }

    ga.setFormat(FLAG_geneFileFormat);
    ga.openGeneFile(FLAG_geneFile.c_str());
    ga.annotate(FLAG_inputFile.c_str(), "test.output.vcf");
#else
    GeneAnnotation ga;
    ga.setFormat(FLAG_geneFileFormat);
    ga.openGeneFile("test.gene.txt");
    ga.openCodonFile("codon.txt");
    ga.openReferenceGenome("test.fa");
    ga.annotate("test.vcf", "test.output.vcf");
#endif
    return 0;
}
