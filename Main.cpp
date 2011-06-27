#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "Argument.h"
#include "IO.h"
#include "StringUtil.h"
#include "LogFile.h"

#include "Codon.h"
#include "GenomeSequence.h"
#include "Gene.h"
#include "GeneFormat.h"
#include "SequenceUtil.h"
#include "FreqTable.h"

void banner(FILE* fp) {
    const char* string =
        "..............................................       \n"
        " ...      G(ene) A(nnotation)                ...     \n"
        "  ...      Xiaowei Zhan, Goncalo Abecasis     ...    \n"
        "   ...      zhanxw@umich.edu                    ...  \n"
        "    ...      Jun 2011                            ... \n"
        "     ................................................\n" 
        "                                                     \n"
        ;
    fputs(string, fp);
};

typedef enum {
    SNP = 0,
    INS,
    DEL,
    MIXED,
    SV,
    UNKNOWN = 99
} VARIATION_TYPE;

typedef enum {
    INSERTION = 0,
    DELETION,
    STRUCTURE_VARIATION,
    UPSTREAM,
    DOWNSTREAM,
    UTR5,
    UTR3,
    INTRON,
    EXON,
    SNV,                    /*SNV contains the following 6 types, it appears when there is no reference.*/
    SYNONYMOUS,
    NONSYNONYMOUS,
    STOP_GAIN,
    STOP_LOSS,
    START_GAIN,
    START_LOSS,
    NORMAL_SPLICE_SITE,
    ESSENTIAL_SPLICE_SITE,
    FRAME_SHIFT,            /* Indel length is not divisible by 3 */
    CODON_GAIN,             /* Insertion length is divisible by 3 */
    CODON_LOSS,             /* Deletion length is divisible by 3 */
    CODON_REGION,           /* Just say the variant is in the Coding Region, used in Structrual Varition*/
    INTROGENIC,
    NONCODING
} AnnotationType;

const char* AnnotationString[]= {
    "Insertion",
    "Deletion",
    "StructuralVariation",
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
    "Frameshift",
    "CodonGain",
    "CodonLoss",
    "CodonRegion",
    "Introgenic",
    "Noncoding"
};


// here we define format or the annotation.
#define FORWARD_STRAND_STRING "(+)"
#define REVERSE_STRAND_STRING "(-)"
#define INFO_SEPARATOR ';'
#define ANNOTATION_START_TAG "ANNO="
#define GENE_SEPARATOR ':'
#define WITHIN_GENE_SEPARATOR '|'
#define WITHIN_GENE_LEFT_DELIM '('
#define WITHIN_GENE_RIGHT_DELIM ')'

/**
 * For each gene, we use AnnotationResult to store all annotation results.
 */
class AnnotationResult{
public:
    void reset() {
        this->data.clear();
        this->freq.clear();
    }
    void clear() {
        this->data.clear();
    };
    std::string toString(){
        // *annotationString += g.name;
        // *annotationString += g.forwardStrand ? FORWARD_STRAND_STRING : REVERSE_STRAND_STRING;
        // annotationString->push_back(WITHIN_GENE_SEPARATOR);
        // *annotationString+=(stringJoin(annotation, WITHIN_GENE_SEPARATOR));

        return stringJoin(this->data, GENE_SEPARATOR);
    };
    void add(const Gene& g){
        if (this->data.size() != 0) {
            // we usually record gene name and its strand first
            fprintf(stderr, "Something weired happen\n");
        }
        this->data.push_back(g.name);
        this->data.push_back(g.forwardStrand ? FORWARD_STRAND_STRING : REVERSE_STRAND_STRING);
    };
    void add(const AnnotationType& t) {
        this->data.push_back(AnnotationString[t]);
        this->freq.add(t);
    };
    // add extra details such as "(CCT/Pro->CAT/His)" to the last element
    template<class T>
    void addDetail(const T& s) {
        unsigned int n = this->data.size();
        assert( n > 0);
        this->data[n-1] += s;
    };
    FreqTable<AnnotationType>& getFreq() {
        return this->freq;
    };
private:
    std::vector<std::string> data;
    FreqTable<AnnotationType> freq;
};
struct GeneAnnotationParam{
    GeneAnnotationParam():
        upstreamRange(50),
        downstreamRange(50),
        spliceIntoExon(3),
        spliceIntoIntron(8) {};
    int upstreamRange;      // upstream def
    int downstreamRange;    // downstream def
    int spliceIntoExon;     // essential splice site def
    int spliceIntoIntron;   // essentail splice site def
};

class GeneAnnotation{
public:
    GeneAnnotation():allowMixedVariation(false){};
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
        fprintf(stdout, "Load gene file %s...\n", geneFileName);
        std::string line;
        std::vector<std::string> fields;
        LineReader lr(geneFileName);
        int totalGene = 0;
        while (lr.readLine(&line) > 0) {
            stringStrip(&line);
            if (line.size()>0 && line[0] == '#' || line.size() == 0) continue; // skip headers and empty lines
            Gene g;
            g.readLine(line.c_str(), this->format);
            this->geneList[g.chr].push_back(g);
            totalGene ++;
        }
        // make sure genes are ordered
        this->sortGene();
        fprintf(stdout, "DONE: %d gene loaded.\n", totalGene);
        LOG << "Gene file " << geneFileName << " loads succeed!\n";
        return;
    };
    void openCodonFile(const char* codonFileName) {
        this->codon.open(codonFileName);
    };
    void openReferenceGenome(const char* referenceGenomeFileName) {
        fprintf(stdout, "Load reference genome %s...\n", referenceGenomeFileName);
        this->gs.open(referenceGenomeFileName);
        fprintf(stdout, "DONE: %d chromosomes and %d bases are loaded.\n", this->gs.size(), this->gs.getGenomeLength());
        LOG << "Reference genome file " << referenceGenomeFileName << " loads succeed!\n";
        return;
    };
    // we take a VCF input file for now
    void annotate(const char* inputFileName, const char* outputFileName){
        // open output file
        FILE* fout = fopen(outputFileName, "wt");
        assert(fout);

        // open input file
        std::vector<std::string> annotationString;
        std::string geneAnnotation;
        LineReader lr(inputFileName);
        std::vector<std::string> field;
        std::string line;
        int totalVariants = 0;
        while(lr.readLine(&line) > 0) {
            if (line.size() == 0 || line[0] == '#') {
                fputs(line.c_str(), fout);
                fputc('\n', fout);
                continue;
            }
            stringTokenize(line, "\t", &field);
            if (field.size() < 4) continue;
            totalVariants++;
            std::string chr = field[0];
            int pos = toInt(field[1]);
            std::string ref = toUpper(field[3]);
            std::string alt = toUpper(field[4]);
            std::vector<unsigned> potentialGeneIdx;
            this->findInRangeGene(field[0], &pos, &potentialGeneIdx);
            // if Introgenic,  we will have (potentialGeneIdx.size() == 0)
            annotationString.clear();
            geneAnnotation.clear();
            for (unsigned int i = 0; i < potentialGeneIdx.size(); i++) {
                this->annotateByGene(potentialGeneIdx[i], chr, pos, ref, alt, &geneAnnotation);
                annotationString.push_back(geneAnnotation);
                // printf("%s %d ref=%s alt=%s has annotation: %s\n",
                //        chr.c_str(), pos,
                //        ref.c_str(), alt.c_str(),
                //        annotationString.c_str());
            }

            // output annotation result
            // VCF info field is the 8th column
            for (unsigned int i = 0; i < field.size(); i++ ){
                if (i) fputc('\t', fout);
                fputs(field[i].c_str(), fout) ;
                if (i == 7) { // 7: the 8th column in 0-based index
                    if (field[i].size() != 0) {
                        fputc(INFO_SEPARATOR, fout);
                    }
                    if (annotationString.size() == 0) {
                        //intergenic
                        fputs(ANNOTATION_START_TAG, fout);
                        fputs(AnnotationString[INTROGENIC], fout);
                    } else {
                        fputs(ANNOTATION_START_TAG, fout);
                        fputs(stringJoin(annotationString, GENE_SEPARATOR).c_str(), fout);
                    }
                }
            }
            fputc('\n', fout);
        }
        // close output
        fclose(fout);
        fprintf(stdout, "DONE: %d varaints are annotated.\n", totalVariants);
        LOG << "Annotate " << inputFileName << " to " << outputFileName << " succeed!\n";
        return;
    };
    void setAnnotationParameter(GeneAnnotationParam& param) {
        this->param = param;
    };
    void printAnnotationFrequency(){
        unsigned int n = this->geneAnnotation.getFreq().size();
        for (unsigned int i = 0; i < n; i++){
            AnnotationType t;
            int count;
            this->geneAnnotation.getFreq().at(i, &t, &count);
            fprintf(stdout, "%s: %d\n", AnnotationString[t], count);
        }
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
        } else if (Codon::isStopCodon(refAAName) && !Codon::isStopCodon(altAAName)) {
            return STOP_LOSS;
        } else if (!Codon::isStopCodon(refAAName) && Codon::isStopCodon(altAAName)) {
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
    /**
     * To annotate for insertion is very similar to annotate SNP, and the only difference is that
     * insertion in the exon could cause frameshift/codon_insertion/codon_deletion
     */
    void annotateIns(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt) {
        Gene& g = this->geneList[chr][geneIdx];
        this->geneAnnotation.add(g);

        // might useful vars.
        int dist2Gene;
        int utrPos, utrLen;
        int exonNum; // which exon
        int codonNum; // which codon
        int codonPos[3] = {0, 0, 0}; // the codon position
        AnnotationType type; // could be one of
        int intronNum; // which intron
        bool isEssentialSpliceSite;

        this->geneAnnotation.add(INSERTION);
        if (g.isUpstream(variantPos, param.upstreamRange, &dist2Gene)) {
            this->geneAnnotation.add(UPSTREAM);
        } else if (g.isDownstream(variantPos, param.upstreamRange, &dist2Gene)) {
            this->geneAnnotation.add(DOWNSTREAM);
        } else if (g.isExon(variantPos, &exonNum)){//, &codonNum, codonPos)) {
            this->geneAnnotation.add(EXON);
            if (!g.isNonCoding()) {
                if (g.is5PrimeUtr(variantPos, &utrPos, &utrLen)) {
                    this->geneAnnotation.add(UTR5);
                } else if (g.is3PrimeUtr(variantPos, &utrPos, &utrLen)) {
                    this->geneAnnotation.add(UTR3);
                } else { // cds part has base change
                    int insertSize = alt.size() - ref.size();
                    if (insertSize % 3 == 0) {
                        this->geneAnnotation.add(CODON_GAIN);
                        std::string s;
                        char triplet[3];
                        std::string aaName;
                        s+= WITHIN_GENE_LEFT_DELIM;
                        for (unsigned int i = ref.size(); i < alt.size();){
                            triplet[0 ] = alt[i++];
                            triplet[1] = alt[i++];
                            triplet[2] = alt[i++];
                            if (!g.forwardStrand)
                                reverseComplementTriplet(triplet);
                            aaName = this->codon.toAA(triplet);
                            s+= aaName;
                        }
                        s+= WITHIN_GENE_RIGHT_DELIM;
                        this->geneAnnotation.addDetail(s);
                    } else {
                        this->geneAnnotation.add(FRAME_SHIFT);
                    }
                }
            } else {
                this->geneAnnotation.add(NONCODING);
            }
            // check splice site
            if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                if (isEssentialSpliceSite)
                    this->geneAnnotation.add(ESSENTIAL_SPLICE_SITE);
                else
                    this->geneAnnotation.add(NORMAL_SPLICE_SITE);
            }
        } else if (g.isIntron(variantPos, &intronNum)) {
            this->geneAnnotation.add(INTRON);
            // check splice site
            if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                if (isEssentialSpliceSite)
                    this->geneAnnotation.add(ESSENTIAL_SPLICE_SITE);
                else
                    this->geneAnnotation.add(NORMAL_SPLICE_SITE);
            }
        } else {
            //annotation.push_back("Intergenic");
        }
    } // end annotateIns(...)
    /**
     * Deletion may across various regions
     * we use std::set to store all regions it came across
     *
     */
    void annotateDel(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt) {
        Gene& g = this->geneList[chr][geneIdx];
        std::set<AnnotationType> annotationSet;

        // preprocessing alt.
        std::string cleanedAlt;
        if (alt == ".") {
            cleanedAlt = "";
        } else {
            cleanedAlt = alt;
        }

        // might useful vars.
        int dist2Gene;
        int utrPos, utrLen;
        int exonNum; // which exon
        int codonNum; // which codon
        int codonPos[3] = {0, 0, 0}; // the codon position
        AnnotationType type; // could be one of
        int intronNum; // which intron
        bool isEssentialSpliceSite;

        // calculate range of the deletion
        // some cases:
        // e.g. ref = AG cleanedAlt = A
        // e.g. ref = AG cleanedAlt = ""
        int delBeg = variantPos + cleanedAlt.size();  // delBeg: inclusive
        int delEnd = variantPos + ref.size(); // delEnd: exclusive

        std::string overlappedCdsBase; // bases in the cds (if any)
        annotationSet.insert(DELETION);
        for (int pos = delBeg; pos < delEnd; pos++) {
            if (g.isUpstream(pos, param.upstreamRange, &dist2Gene)) {
                annotationSet.insert(UPSTREAM);
            } else if (g.isDownstream(pos, param.upstreamRange, &dist2Gene)) {
                annotationSet.insert(DOWNSTREAM);
            } else if (g.isExon(pos, &exonNum)){//, &codonNum, codonPos)) {
                annotationSet.insert(EXON);
                if (!g.isNonCoding()) {
                    if (g.is5PrimeUtr(pos, &utrPos, &utrLen)) {
                        annotationSet.insert(UTR5);
                    } else if (g.is3PrimeUtr(pos, &utrPos, &utrLen)) {
                        annotationSet.insert(UTR3);
                    } else { // cds part has base change
                        overlappedCdsBase.push_back(ref[ pos - delBeg ]);
                    }
                } else {
                    annotationSet.insert(NONCODING);
                }
                // check splice site
                if (g.isSpliceSite(pos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                    if (isEssentialSpliceSite)
                        annotationSet.insert(ESSENTIAL_SPLICE_SITE);
                    else
                        annotationSet.insert(NORMAL_SPLICE_SITE);
                }
            } else if (g.isIntron(pos, &intronNum)) {
                annotationSet.insert(INTRON);
                // check splice site
                if (g.isSpliceSite(pos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                    if (isEssentialSpliceSite)
                        annotationSet.insert(ESSENTIAL_SPLICE_SITE);
                    else
                        annotationSet.insert(NORMAL_SPLICE_SITE);
                }
            } else {
                //annotation.push_back("Intergenic");
            }
        } // end for
        // check how many codon in cds are delete
        if (overlappedCdsBase.size() > 0) {
            if (overlappedCdsBase.size() % 3 == 0) {
                annotationSet.insert(CODON_LOSS);
            } else {
                annotationSet.insert(FRAME_SHIFT);
            }
        }
        // store all existing annotation
        this->geneAnnotation.add(g);
        for (std::set<AnnotationType>::const_iterator it = annotationSet.begin();
             it != annotationSet.end();
             it++) {
            this->geneAnnotation.add(*it);
        };
    };
    /**
     * SV is the most complex scenario. fully support this is an ongoing work.
     * We will just annotation the region in the rough scale.
     */
    void annotateSV(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt) {
        Gene& g = this->geneList[chr][geneIdx];
        this->geneAnnotation.add(g);

        // might useful vars.
        int dist2Gene;
        int utrPos, utrLen;
        int exonNum; // which exon
        int codonNum; // which codon
        int codonPos[3] = {0, 0, 0}; // the codon position
        AnnotationType type; // could be one of
        int intronNum; // which intron
        bool isEssentialSpliceSite;
        this->geneAnnotation.add(STRUCTURE_VARIATION);
        if (g.isUpstream(variantPos, param.upstreamRange, &dist2Gene)) {
            this->geneAnnotation.add(UPSTREAM);
        } else if (g.isDownstream(variantPos, param.upstreamRange, &dist2Gene)) {
            this->geneAnnotation.add(DOWNSTREAM);
        } else if (g.isExon(variantPos, &exonNum)){//, &codonNum, codonPos)) {
            this->geneAnnotation.add(EXON);
            if (!g.isNonCoding()) {
                if (g.is5PrimeUtr(variantPos, &utrPos, &utrLen)) {
                    this->geneAnnotation.add(UTR5);
                } else if (g.is3PrimeUtr(variantPos, &utrPos, &utrLen)) {
                    this->geneAnnotation.add(UTR3);
                } else { // cds part has base change
                    this->geneAnnotation.add(CODON_REGION);
                }
            } else{
                this->geneAnnotation.add(NONCODING);
            }
            // check splice site
            if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                if (isEssentialSpliceSite)
                    this->geneAnnotation.add(ESSENTIAL_SPLICE_SITE);
                else
                    this->geneAnnotation.add(NORMAL_SPLICE_SITE);
            }
        } else if (g.isIntron(variantPos, &intronNum)) {
            this->geneAnnotation.add(INTRON);
            // check splice site
            if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                if (isEssentialSpliceSite)
                    this->geneAnnotation.add(ESSENTIAL_SPLICE_SITE);
                else
                    this->geneAnnotation.add(NORMAL_SPLICE_SITE);
            }
        } else {
            //annotation.push_back("Intergenic");
        }
    }
    void annotateSNP(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt) {
        Gene& g = this->geneList[chr][geneIdx];
        this->geneAnnotation.add(g);

        // might useful vars.
        int dist2Gene;
        int utrPos, utrLen;
        int exonNum; // which exon
        int codonNum; // which codon
        int codonPos[3] = {0, 0, 0}; // the codon position
        AnnotationType type; // could be one of
        int intronNum; // which intron
        bool isEssentialSpliceSite;

        if (g.isUpstream(variantPos, param.upstreamRange, &dist2Gene)) {
            this->geneAnnotation.add(UPSTREAM);
        } else if (g.isDownstream(variantPos, param.upstreamRange, &dist2Gene)) {
            this->geneAnnotation.add(DOWNSTREAM);
        } else if (g.isExon(variantPos, &exonNum)){//, &codonNum, codonPos)) {
            this->geneAnnotation.add(EXON);
            if (!g.isNonCoding()) {
                if (g.is5PrimeUtr(variantPos, &utrPos, &utrLen)) {
                    this->geneAnnotation.add(UTR5);
                } else if (g.is3PrimeUtr(variantPos, &utrPos, &utrLen)) {
                    this->geneAnnotation.add(UTR3);
                } else { // cds part has base change
                    if (g.calculateCodonPosition(variantPos, &codonNum, codonPos)) {
                        char refTriplet[3];
                        char altTriplet[3];
                        std::string refAAName;
                        std::string altAAName;
                        AnnotationType annotationType;
                        // when reference genome is provided
                        if (this->gs.size() >0){
                            this->fillTriplet(chr, variantPos, codonPos, g.forwardStrand, ref, alt, refTriplet, altTriplet);

                            refAAName = this->codon.toAA(refTriplet);
                            altAAName = this->codon.toAA(altTriplet);
                            annotationType = this->determineSNVType(refAAName, altAAName, codonNum);

                            this->geneAnnotation.add(annotationType);
                            std::string s;
                            s += WITHIN_GENE_LEFT_DELIM;
                            s += refTriplet[0];
                            s += refTriplet[1];
                            s += refTriplet[2];
                            s += "/";
                            s += refAAName;
                            s += "->";
                            s += altTriplet[0];
                            s += altTriplet[1];
                            s += altTriplet[2];
                            s += "/";
                            s += altAAName;
                            s += WITHIN_GENE_RIGHT_DELIM;
                            this->geneAnnotation.addDetail(s);
                        } else {
                            this->geneAnnotation.add(SNV);
                        }
                    } else {
                        this->geneAnnotation.add(SNV);
                    }
                }
            } else{
                this->geneAnnotation.add(NONCODING);
            }
            // check splice site
            if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                if (isEssentialSpliceSite)
                    this->geneAnnotation.add(ESSENTIAL_SPLICE_SITE);
                else
                    this->geneAnnotation.add(NORMAL_SPLICE_SITE);
            }
        } else if (g.isIntron(variantPos, &intronNum)) {
            this->geneAnnotation.add(INTRON);
            // check splice site
            if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
                if (isEssentialSpliceSite)
                    this->geneAnnotation.add(ESSENTIAL_SPLICE_SITE);
                else
                    this->geneAnnotation.add(NORMAL_SPLICE_SITE);
            }
        } else {
            //annotation.push_back("Intergenic");
        }
    }
    /**
     *
     */
    void annotateByGene(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt,
                        std::string* annotationString){
        assert(annotationString);
        annotationString->clear();
        this->geneAnnotation.clear();
        // check VARATION_TYPE
        switch(determineVariationType(ref, alt)) {
        case SNP:
            this->annotateSNP(geneIdx, chr, variantPos, ref, alt);
            break;
        case INS:
            this->annotateIns(geneIdx, chr, variantPos, ref, alt);
            break;
        case DEL:
            this->annotateDel(geneIdx, chr, variantPos, ref, alt);
            break;
        case MIXED:
            if (this->allowMixedVariation) {
                // annotate multiple gene
                // to finished.
            } else {
                // only annotate the first variation
                std::string singleAlt = alt;
                int commaPos = alt.find(',');
                singleAlt = alt.substr(0, commaPos);
                this->annotateByGene(geneIdx, chr, variantPos, ref, singleAlt, annotationString);
            };
            break;
        case SV:
            this->annotateSV(geneIdx, chr, variantPos, ref, alt);
            break;
        case UNKNOWN:
        default:
            break;
        };

        annotationString->assign(this->geneAnnotation.toString());
        return;
    };
    /**
     * @return the variation type depending on the first entry in the alt field
     */
    VARIATION_TYPE determineVariationType(const std::string& ref, const std::string& alt) {
        if (alt.find(',') != std::string::npos) {
            return MIXED;
        }
        // NOTE: a single "." is used for deletion; but "G." can represent uncertain breakpoint
        //       we will check the former case first.
        if (alt == ".") {
            return DEL;
        }

        const char* ALLOWED_BASE = "ACGT";
        unsigned int refLen = ref.size();
        unsigned int altLen = alt.size();
        if (alt.find_first_not_of(ALLOWED_BASE) != std::string::npos) {
            // NOTE: SV usually contain "[" or "]" for rearrangment
            //       ">", "<" for haplotypes or large deletion/insertion.
            return SV;
        }
        if (refLen == altLen) {
            if (refLen == 1){
                return SNP;
            } else {
                return UNKNOWN;
            }
        } else if (refLen > altLen) {
            return DEL;
        } else if (refLen < altLen) {
            return INS;
        }
        return UNKNOWN;
    };
private:
    GeneAnnotationParam param;
    std::map <std::string, std::vector<Gene> > geneList;
    std::string annotation;
    GenomeSequence gs;
    Codon codon;
    GeneFormat format;
    bool allowMixedVariation;       // VCF ALT field may have more than one variation e..g A,C

    AnnotationResult geneAnnotation;
    // FreqTable<std::string> annotationTypeFreq; // annotation type frequency
    // FreqTable<std::string> baseFreq;           // base change frequency
    // FreqTable<std::string> codonFreq;          // codon change frequency
    // FreqTable<int> indelLengthFreq; // for insertion, the value is positive; for deletion, positive
};

int main(int argc, char *argv[])
{
    banner(stdout);

    BEGIN_PARAMETER_LIST(pl)
        ADD_STRING_PARAMETER(pl, inputFile, "-i", "Specify input VCF file")
        ADD_STRING_PARAMETER(pl, outputFile, "-o", "Specify output VCF file")
        ADD_STRING_PARAMETER(pl, geneFile, "-g", "Specify gene file")
        ADD_STRING_PARAMETER(pl, referenceFile, "-r", "Specify reference genome position")
        ADD_STRING_PARAMETER(pl, geneFileFormat, "-f", "Specify gene file format (default: refFlat, other options knownGene)")
        END_PARAMETER_LIST(pl)
        ;

    pl.Read(argc, argv);
    if (FLAG_geneFileFormat.size() == 0) {
        FLAG_geneFileFormat = "refFlat";
    }

    if (FLAG_inputFile.size() == 0) {
        pl.Help();
        fprintf(stderr, "Please specify input file\n");
        exit(1);
    }
    if (FLAG_outputFile.size() == 0) {
        pl.Help();
        fprintf(stderr, "Please specify output file\n");
        exit(1);
    }
    if (FLAG_geneFile.size() == 0) {
        pl.Help();
        fprintf(stderr, "Please specify gene file\n");
        exit(1);
    }

#if 1
    std::string logFileName = FLAG_outputFile + ".log";
    LOG_START(logFileName.c_str());
    LOG_START_TIME;
    LOG_PARAMETER(pl);
    GeneAnnotation ga;
    pl.Status();
    if (FLAG_referenceFile.size() != 0) {
        ga.openReferenceGenome(FLAG_referenceFile.c_str());
        ga.openCodonFile("codon.txt");
    }

    ga.setFormat(FLAG_geneFileFormat);
    ga.openGeneFile(FLAG_geneFile.c_str());
    ga.annotate(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
    ga.printAnnotationFrequency();

    LOG_END_TIME;
    LOG_END ;
#else
    // debug purpose
    GeneAnnotation ga;
    ga.setFormat(FLAG_geneFileFormat);
    ga.openGeneFile("test.gene.txt");
    ga.openCodonFile("codon.txt");
    ga.openReferenceGenome("test.fa");
    ga.annotate("test.vcf", "test.output.vcf");
#endif
    printf("Annotation succeed!\n");
    return 0;
}
