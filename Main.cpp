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
#include <algorithm>

#include "Argument.h"
#include "IO.h"
#include "StringUtil.h"

#include "Codon.h"
#include "GenomeSequence.h"
#include "Gene.h"
#include "GeneFormat.h"
#include "SequenceUtil.h"

typedef enum {
    SNP = 0,
    INS,
    DEL,
    MIXED,
    SV,
    UNKNOWN = 99
} VARIATION_TYPE;

typedef enum {
    UPSTREAM = 0,
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
    FRAME_SHIFT,
    CODON_GAIN,
    CODON_LOSS,
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
    "Frameshift",
    "CodonGain",
    "CodonLoss",
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
        // open output file
        FILE* fout = fopen(outputFileName, "wt");
        assert(fout);
        
        // open input file
        std::string annotationString;
        LineReader lr(inputFileName);
        std::vector<std::string> field;
        std::string line;
        while(lr.readLine(&line) > 0) {
            if (line.size() == 0 || line[0] == '#') {
                fputs(line.c_str(), fout);
                fputc('\n', fout);
                continue;
            }
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
                    if (annotationString.size() == 0) {
                        //intergenic
                        fputs(";ANNO=", fout);
                        fputs(AnnotationString[INTROGENIC], fout);
                    } else {
                        fputs(";ANNO=", fout);
                        fputs(annotationString.c_str(), fout);
                    }
                }
            }
            fputc('\n', fout);
        }
        // close output
        fclose(fout);

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
    void annotateIns(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt,
                     std::string* annotationString) {
    }
    void annotateDel(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt,
                     std::string* annotationString) {
    }
    void annotateSV(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt,
                     std::string* annotationString) {
    }
    void annotateSNP(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt,
                     std::string* annotationString) {
        Gene& g = this->geneList[chr][geneIdx];
        std::vector<std::string> annotation;

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
            annotation.push_back(AnnotationString[UPSTREAM]);
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
                            this->fillTriplet(chr, variantPos, codonPos, g.forwardStrand, ref, alt, refTriplet, altTriplet);

                            refAAName = this->codon.toAA(refTriplet);
                            altAAName = this->codon.toAA(altTriplet);
                            annotationType = this->determineSNVType(refAAName, altAAName, codonNum);

                            std::string s;
                            s = AnnotationString[annotationType];
                            s += '(';
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
                            s += ')';
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
        *annotationString += g.name;
        *annotationString += g.forwardStrand ? "(+)" : "(-)"; 
        annotationString->push_back('|');
        *annotationString+=(stringJoin(annotation, "|"));
    }
    /**
     * 
     */
    void annotateByGene(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt,
                        std::string* annotationString){
        // check VARATION_TYPE
        switch(determineVariationType(ref, alt)) {
        case SNP:
            this->annotateSNP(geneIdx, chr, variantPos, ref, alt, annotationString);
            break;
        case INS:
            this->annotateIns(geneIdx, chr, variantPos, ref, alt, annotationString);
            break;
        case DEL:
            this->annotateDel(geneIdx, chr, variantPos, ref, alt, annotationString);
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
            this->annotateSV(geneIdx, chr, variantPos, ref, alt, annotationString);
            break;
        case UNKNOWN:
        default:
            break;
        };
        return;
    };
    /**
     * @return the variation type depending on the first entry in the alt field
     */
    VARIATION_TYPE determineVariationType(const std::string& ref, const std::string& alt) {
        if (alt.find(',') != std::string::npos) {
            return MIXED;
        }
        const char* BASE = "ACGT";
        unsigned int refLen = ref.size();
        unsigned int altLen = alt.size();
        if (alt.find_first_not_of(BASE) != std::string::npos) {
            //contain 
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
    bool allowMixedVariation; // ALT may has more than one variation e..g A,C
};
int main(int argc, char *argv[])
{
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
    GeneAnnotation ga;
    pl.Status();
    if (FLAG_referenceFile.size() != 0) {
        ga.openReferenceGenome(FLAG_referenceFile.c_str());
        ga.openCodonFile("codon.txt");
    }

    ga.setFormat(FLAG_geneFileFormat);
    ga.openGeneFile(FLAG_geneFile.c_str());
    ga.annotate(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
#else
    // debug purpose
    GeneAnnotation ga;
    ga.setFormat(FLAG_geneFileFormat);
    ga.openGeneFile("test.gene.txt");
    ga.openCodonFile("codon.txt");
    ga.openReferenceGenome("test.fa");
    ga.annotate("test.vcf", "test.output.vcf");
#endif
    return 0;
}
