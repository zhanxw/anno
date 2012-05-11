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
#include "StringTemplate.h"

void banner(FILE* fp) {
  const char* string =
      "..............................................       \n"
      " ...      Anno(tation)                       ...     \n"
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
  INTERGENIC,
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
  "Intergenic",
  "Noncoding"
};


// here we define format or the annotation.
#define FORWARD_STRAND_STRING "+"
#define REVERSE_STRAND_STRING "-"
#define WITHIN_GENE_LEFT_DELIM '('
#define WITHIN_GENE_RIGHT_DELIM ')'
#define ANNOTATION_START_TAG "ANNO="
#define ANNOTATION_START_TAG_FULL "ANNOFULL="
#define VCF_INFO_SEPARATOR ';'
#define WITHIN_GENE_SEPARATOR ':'
// #define GENE_SEPARATOR '|'



/**
 * priority are given by lineNo
 * small number -> high priority (more important)
 * contain a relationship between AnnotationType and int(priority)
 */
struct Priority{
 public:
  int open(const char* fileName) {
    // Load priority.txt
    this->priorityIdx = 0;
    this->priorityInt2Str.clear();
    this->priorityStr2Int.clear();
    LineReader lr(fileName);
    std::vector<std::string> fd;

    while (lr.readLineBySep(&fd, " \t")){
      if (fd.size() == 0) continue;
      if (fd[0][0] == '#') continue;
      if (fd[0].size() == 0) continue;
      priorityIdx ++;
      // fprintf(stderr, "add priority [%s]\n", fd[0].c_str());
      priorityInt2Str[priorityIdx] = fd[0];
      priorityStr2Int[fd[0]] = priorityIdx;
    }
    return priorityIdx;
  };
  int getPriority(const AnnotationType& t) const{
    std::map<std::string, int>::const_iterator it;
    it = this->priorityStr2Int.find( AnnotationString[t] );
    if (it == this->priorityStr2Int.end()) {
      fprintf(stderr, "Cannot find annotation type %s from priority files!", AnnotationString[t]);
      return -1;
    } else {
      return it->second;
    }
  };

  std::string getAnnotationString (const int& i) const{
    std::map<int, std::string>::const_iterator it;
    it = this->priorityInt2Str.find( i );
    if (it == this->priorityInt2Str.end()) {
      fprintf(stderr, "Cannot find priority %d from priority files!", i);
      return "";
    } else {
      return it->second;
    }
  };
  /**
   * @return  @param priority
   */
  std::string toString(const int priority) const {
    std::string s;
    std::map<int, std::string>::const_iterator it;
    it = this->priorityInt2Str.find(priority);
    if (it == this->priorityInt2Str.end()) {
      return s;
    } else{
      return it->second;
    }
  };

  static int getLeastPriority(){
    return 9999;
  }
 private:
  int priorityIdx;
  std::map<int, std::string> priorityInt2Str;
  std::map<std::string, int> priorityStr2Int;
}; // end Priority

/**
 * For each gene, we use AnnotationResult to store all annotation results.
 */
class AnnotationResult{
 public:
  void clear() {
    this->geneName.clear();
    this->type.clear();
    this->detail.clear();
    this->topPriorityIndex = -1;
  };
#if 0
  /**
   * will need to refactor using annotation output
   */
  std::string toString() const{
    fprintf(stderr, "Deprecated!\n");
    std::string s;
    s += this->geneName;
    s += WITHIN_GENE_SEPARATOR;
    s += isForwardStrand ? FORWARD_STRAND_STRING : REVERSE_STRAND_STRING;

    for (int i = 0; i < type.size(); i++) {
      s += WITHIN_GENE_SEPARATOR;
      s += AnnotationString[type[i]];
      // for synonymous or non-synonymous, we may output details
      std::map<AnnotationType, std::string>::const_iterator iter;
      iter = this->detail.find(type[i]);
      if (iter != this->detail.end()) {
        s += iter->second;
      }
    }
    return s;
  };
#endif
  std::string toString(StringTemplate* t) const{
    std::vector<std::string> type;
    std::string s;
    for (int i = 0; i < this->type.size(); i++) {
      s =  AnnotationString[this->type[i]];
      // for synonymous or non-synonymous, we may output details
      std::map<AnnotationType, std::string>::const_iterator iter;
      iter = this->detail.find(this->type[i]);
      if (iter != this->detail.end()) {
        s += iter->second;
      }
      type.push_back(s);
    }
    t->add("GENE_NAME", this->geneName);
    t->add("GENE_STRAND", isForwardStrand ? FORWARD_STRAND_STRING : REVERSE_STRAND_STRING);
    t->add("TYPE", type);
    
    t->translate(&s);
    return s;
  };

  void add(const Gene& g){
    if (this->type.size() != 0) {
      // we usually record gene name and its strand first
      fprintf(stderr, "Something weired happen\n");
    }
    this->geneName = g.name;
    this->isForwardStrand = g.forwardStrand;
    //this->data.push_back(g.forwardStrand ? FORWARD_STRAND_STRING : REVERSE_STRAND_STRING);
  };
  void add(const AnnotationType& t) {
    this->type.push_back(t);
  };
  // add extra details such as "(CCT/Pro->CAT/His)" to the last element
  void addDetail(const AnnotationType& t, const std::string& s) {
    this->detail[t] = s;
  };

  //////////////////////////////////////////////////////////////////////
  // priority related functions
  /**
   * set @param topPriorityIndex to the index in this->type with top priority
   */
  int calculateTopPriority(const Priority& p) {
    if (this->topPriorityIndex >= 0)
      return this->topPriorityIndex;

    this->topPriorityIndex = -1;
    int highest = Priority::getLeastPriority();
    for (unsigned int i = 0; i != type.size(); i++) {
      int t = p.getPriority(type[i]);
      if (t <= highest){
        highest = t;
        this->topPriorityIndex = i;
      }
    }
    return this->topPriorityIndex;
  };

  // /**
  //  * @return top priority string.
  //  * NOTE: should call calculateTopPriority first!
  //  */
  // std::string getTopPriorityString() const{
  //   assert( this->topPriorityIndex >= 0);
  //   std::string s = AnnotationString[this->topPriorityIndex]; // #p.toString(highest);
  //   s += WITHIN_GENE_SEPARATOR;
  //   s += this->geneName;
  //   return s;
  // };

  AnnotationType getTopPriorityType() const{
    assert( this->topPriorityIndex >= 0);
    return this->type[this->topPriorityIndex];
  };

  //////////////////////////////////////////////////////////////////////
  // getters
  const std::string& getGeneName() const{
    return this->geneName;
  };
  
  const std::vector<AnnotationType>& getType() const{
    return this->type;
  };


 private:
  std::string geneName;
  bool isForwardStrand;
  std::vector<AnnotationType> type;
  std::map<AnnotationType, std::string> detail;
  int topPriorityIndex;  // this->type[this->topPriorityIndex] has top priority; <0, means unknown
}; // end AnnotationResult

class AnnotationOutput{
 public:
  AnnotationOutput():
      topPriorityTemplate("$(TOP_TYPE):$[$(ALL_TOP_GENE) |]"),
      geneTemplate("$(GENE_NAME):$(GENE_STRAND):$[$(TYPE) :]"),
      fullTemplate("$[$(GENE_ANNOTATION) |]"),
      annotationResult(NULL),
      priority(NULL)
  {
  };
  void setAnnotationResult(const std::vector<AnnotationResult>& r) {
    this->annotationResult = &r;
    this->buildKeywordDict();
  };
  void setPriority(const Priority& p) {
    this->priority = &p;
  };
  void buildKeywordDict() {
    const std::vector<AnnotationResult>& r = *this->annotationResult;
    std::vector<std::string> res;
    int highestPriority = Priority::getLeastPriority();
    unsigned int highestIdx = 0;
    for (unsigned int i = 0; i != r.size(); i++) {
      int p = this->priority->getPriority(r[i].getTopPriorityType());
      if (p < highestPriority) {
        highestPriority = p;
        highestIdx = i;
        res.clear();
        res.push_back(r[i].getGeneName());
      } else if (p == highestPriority) {
        res.push_back(r[i].getGeneName());
      }
    }
    topPriorityTemplate.add("TOP_GENE", res.size() ? res[0]: AnnotationString[INTERGENIC]);
    topPriorityTemplate.add("ALL_TOP_GENE", res);
    topPriorityTemplate.add("TOP_TYPE", this->priority->toString(highestPriority).c_str());

    std::vector<std::string> geneTemplate;
    std::string perGeneAnnotation;
    for (unsigned int i = 0; i != r.size(); i++){
      perGeneAnnotation = r[i].toString(&this->geneTemplate);
      geneTemplate.push_back(perGeneAnnotation);
    }
    fullTemplate.add("GENE_ANNOTATION", geneTemplate);
  };
  // Format:
  //  Most_priority:gene1|gene2
  std::string getTopPriorityAnnotation(const std::vector<AnnotationResult>& r, const Priority& p) const{
    if (r.size() == 0) {
      return AnnotationString[INTERGENIC];
    };
    std::string s;
    if (this->topPriorityTemplate.translate(&s)) {
      fprintf(stderr, "topPriorityTemplate failed translation!\n");
    }
    return s;
  };

  std::string getFullAnnotation(const std::vector<AnnotationResult>& r){
    if (r.size() == 0) {
      return AnnotationString[INTERGENIC];
    };
    std::string s;
    if (this->fullTemplate.translate(&s)) {
      fprintf(stderr, "fullTemplate failed translation!\n");
    }
    return s;
  };
 private:
  StringTemplate topPriorityTemplate;
  StringTemplate geneTemplate;
  StringTemplate fullTemplate;
  const std::vector<AnnotationResult>* annotationResult;
  const Priority* priority;
}; // end AnnotationOutput

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
    fprintf(stdout, "Load codon file %s...\n", codonFileName);
    this->codon.open(codonFileName);
    fprintf(stdout, "DONE: codon file loaded.\n");
    LOG << "Codon file " << codonFileName << " loads succeed!\n";
    return;
  };
  void openReferenceGenome(const char* referenceGenomeFileName) {
    fprintf(stdout, "Load reference genome %s...\n", referenceGenomeFileName);
    this->gs.open(referenceGenomeFileName);
    fprintf(stdout, "DONE: %d chromosomes and %ld bases are loaded.\n", this->gs.size(), this->gs.getGenomeLength());
    LOG << "Reference genome file " << referenceGenomeFileName << " loads succeed!\n";
    return;
  };
  void openPriorityFile(const char* fileName) {
    fprintf(stdout, "Load priority file %s...\n", fileName);
    int ret = this->priority.open(fileName);
    fprintf(stderr, "DONE: %d priority annotation types loaded.\n", ret);
    LOG << "Priority file " << fileName << " load succeed!\n";

    this->outputter.setPriority(this->priority);
    return;
  };
  /**
   * parse input file, and call annotate, and then output results
   */
  void annotateVCF(const char* inputFileName, const char* outputFileName){
    // open output file
    FILE* fout = fopen(outputFileName, "wt");
    assert(fout);

    // open input file
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
      if (field.size() < 5) continue;
      totalVariants++;
      std::string chr = field[0];
      int pos = toInt(field[1]);
      std::string ref = toUpper(field[3]);
      std::string alt = toUpper(field[4]);

      // real part of annotation
      annotate(chr, pos, ref, alt, &annotationResult);
      outputter.setAnnotationResult(annotationResult);
      
      // output annotation result
      // VCF info field is the 8th column
      for (unsigned int i = 0; i < field.size(); i++ ){
        if (i) fputc('\t', fout);
        fputs(field[i].c_str(), fout) ;
        if (i == 7) { // 7: the 8th column in 0-based index
          if (field[i].size() != 0) {
            fputc(VCF_INFO_SEPARATOR, fout);
          }
          fputs(ANNOTATION_START_TAG, fout);
          fputs(outputter.getTopPriorityAnnotation(annotationResult, this->priority).c_str(), fout);
          fputc(VCF_INFO_SEPARATOR, fout);
          fputs(ANNOTATION_START_TAG_FULL, fout);
          fputs(outputter.getFullAnnotation(annotationResult).c_str(), fout);
        }
      }
      fputc('\n', fout);
    }
    // close output
    fclose(fout);
    fprintf(stdout, "DONE: %d varaints are annotated.\n", totalVariants);
    LOG << "Annotate " << inputFileName << " to " << outputFileName << " succeed!\n";

    // output stats
    this->outputAnnotationStats(outputFileName);
  }
  void annotatePlain(const char* inputFileName, const char* outputFileName){
    // open output file
    FILE* fout = fopen(outputFileName, "wt");
    assert(fout);

    // open input file
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
      std::string ref = toUpper(field[2]);
      std::string alt = toUpper(field[3]);

      // real part of annotation
      annotate(chr, pos, ref, alt, &annotationResult);
      outputter.setAnnotationResult(annotationResult);

      // output annotation result
      // VCF info field is the 8th column
      for (unsigned int i = 0; i < field.size(); i++ ){
        fputs(field[i].c_str(), fout) ;
        fputs("\t", fout);
      }
      fputs(outputter.getTopPriorityAnnotation(annotationResult, this->priority).c_str(), fout);
      fputs("\t", fout);
      fputs(outputter.getFullAnnotation(annotationResult).c_str(), fout);
      fputc('\n', fout);
    }
    // close output
    fclose(fout);
    fprintf(stdout, "DONE: %d varaints are annotated.\n", totalVariants);
    LOG << "Annotate " << inputFileName << " to " << outputFileName << " succeed!\n";

    // output stats
    this->outputAnnotationStats(outputFileName);
  }
  // plink associated results (.assoc files)
  // NOTE: ref and alt allelel need to get from reference
  //    1   1:196404269  196404269    A       19        0    G         6.61      0.01014           NA
  void annotatePlink(const char* inputFileName, const char* outputFileName){
    // open output file
    FILE* fout = fopen(outputFileName, "wt");
    assert(fout);

    // open input file
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
      stringNaturalTokenize(line, " ", &field);
      if (field.size() < 10) continue;
      totalVariants++;
      std::string chr = field[0];
      int pos = toInt(field[2]);
      std::string ref = toUpper(field[3]);
      std::string alt = toUpper(field[6]);

      // determine ref base from reference
      bool refMatchRef = true;
      for (int i = 0; i < ref.size(); i++ ) {
        if (ref[i] != gs[chr][pos - 1 + i]) {
          refMatchRef = false;
          break;
        }
      }
      if (!refMatchRef) {
        bool altMatchRef = true;
        for (int i = 0; i < alt.size(); i++ ) {
          if (alt[i] != gs[chr][pos - 1 + i]) {
            altMatchRef = false;
            break;
          }
        }
        if (!altMatchRef) {
          fprintf(stderr, "Ref [%s] and alt [%s] does not match reference: %s:%d\n", ref.c_str(), alt.c_str(), chr.c_str(), pos);
          continue;
        } else {
          std::swap(ref, alt);
        }
      }

      // real part of annotation
      annotate(chr, pos, ref, alt, &annotationResult);
      outputter.setAnnotationResult(annotationResult);

      // output annotation result
      // VCF info field is the 8th column
      for (unsigned int i = 0; i < field.size(); i++ ){
        fputs(field[i].c_str(), fout) ;
        fputs("\t", fout);
      }
      fputs(outputter.getTopPriorityAnnotation(annotationResult, this->priority).c_str(), fout);
      fputs("\t", fout);
      fputs(outputter.getFullAnnotation(annotationResult).c_str(), fout);
      fputc('\n', fout);
    }
    // close output
    fclose(fout);
    fprintf(stdout, "DONE: %d varaints are annotated.\n", totalVariants);
    LOG << "Annotate " << inputFileName << " to " << outputFileName << " succeed!\n";

    // output stats
    this->outputAnnotationStats(outputFileName);
  }
  void outputAnnotationStats(const char* outputFileName) {
    // output frequency files
    std::string fn = outputFileName;


    // output annotation frequency (all types of annotation)
    std::string ofs = fn+".anno.frq";
    this->printAnnotationFrequency(ofs.c_str());
    fprintf(stdout, "DONE: Generated frequency of each annotype type in [ %s ].\n", ofs.c_str());
    LOG << "Generate frequency of each annotation type in " << ofs << " succeed!\n";

    // output annotation frequency
    ofs = fn+".top.anno.frq";
    this->printTopPriorityAnnotationFrequency(ofs.c_str());
    fprintf(stdout, "DONE: Generated frequency of each highest priority annotation type in [ %s ].\n", ofs.c_str());
    LOG << "Generate frequency of high priority for highest priority annotation type in " << ofs << " succeed!\n";
    
    // output base change frequency
    ofs = fn+".base.frq";
    this->printBaseChangeFrequency(ofs.c_str());
    fprintf(stdout, "DONE: Generated frequency of each base change in [ %s ].\n", ofs.c_str());
    LOG << "Generate frequency of each base change in " << ofs << " succeed!\n";

    // output codon change frequency
    ofs = fn+".codon.frq";
    this->printCodonChangeFrequency(ofs.c_str());
    fprintf(stdout, "DONE: Generated frequency of each codon change in [ %s ].\n", ofs.c_str());
    LOG << "Generate frequency of each codon change in " << ofs << " succeed!\n";

    // output indel length frequency
    ofs = fn+".indel.frq";
    this->printIndelLengthFrequency(ofs.c_str());
    fprintf(stdout, "DONE: Generated frequency of indel length in [ %s ].\n", ofs.c_str());
    LOG << "Generate frequency of indel length in " << ofs << " succeed!\n";
  };
  // annotation will find all overlapping gene, and call annotateByGene for each gene
  void annotate(const std::string& chrom,
                const int pos,
                const std::string& ref,
                const std::string& altParam,
                std::vector<AnnotationResult>* annotationResult) {
    // check VARATION_TYPE
    std::string alt = altParam;
    VARIATION_TYPE type = determineVariationType(ref, alt);
    if (type == MIXED) {
      // only annotate the first variation
      int commaPos = alt.find(',');
      alt = altParam.substr(0, commaPos);
      type = determineVariationType(ref, alt);
    }

    // find near target genes
    std::vector<unsigned> potentialGeneIdx;
    this->findInRangeGene(chrom, pos, &potentialGeneIdx);
    // if Intergenic,  we will have (potentialGeneIdx.size() == 0)
    this->annotationResult.clear();
    AnnotationResult annotationPerGene;

    // determine variation type
    for (unsigned int i = 0; i < potentialGeneIdx.size(); i++) {
      annotationPerGene.clear();
      this->annotateByGene(potentialGeneIdx[i], chrom, pos, ref, alt, type, &annotationPerGene);
      this->annotationResult.push_back(annotationPerGene);
    }

    // record frquency info
    if (!this->annotationResult.size()) { // intergenic
      this->annotationTypeFreq.add(INTERGENIC);
      this->topPriorityAnnotationTypeFreq.add(INTERGENIC);
    } else {
      int topPriorityIndex = -1;
      int topPriority = Priority::getLeastPriority();
      for (unsigned int i = 0; i != this->annotationResult.size(); i++) {
        // calculate high priority for each gene
        this->annotationResult[i].calculateTopPriority(this->priority);

        //
        const AnnotationType& t = this->annotationResult[i].getTopPriorityType();
        int p = this->priority.getPriority(t);
        if ( p < topPriority) {
          topPriorityIndex = i;
          topPriority = p;
        }
        for (unsigned int j = 0; j < this->annotationResult[i].getType().size(); j++ ){
          this->annotationTypeFreq.add( this->annotationResult[i].getType() [j] );
        }
      }
      this->topPriorityAnnotationTypeFreq.add(this->annotationResult[topPriorityIndex].getTopPriorityType());
    }
    switch (type) {
      case SNP:
        this->baseFreq.add(ref + "->" +alt);
        break;
      case INS:
      case DEL:
        this->indelLengthFreq.add(calculateIndelLength(ref, alt));
      default:
        break;
    }
    return;
  };
  void setAnnotationParameter(GeneAnnotationParam& param) {
    this->param = param;
  };
  void printAnnotationFrequency(const char* fileName){
    FILE* fp = fopen(fileName, "wt");
    assert(fp);
    unsigned int n = this->annotationTypeFreq.size();
    for (unsigned int i = 0; i < n; i++){
      AnnotationType t;
      int freq;
      this->annotationTypeFreq.at(i, &t, &freq);
      fprintf(fp, "%s\t%d\n", AnnotationString[t], freq);
    }
    fclose(fp);
  };
  void printTopPriorityAnnotationFrequency(const char* fileName){
    FILE* fp = fopen(fileName, "wt");
    assert(fp);
    unsigned int n = this->topPriorityAnnotationTypeFreq.size();
    for (unsigned int i = 0; i < n; i++){
      AnnotationType t;
      int freq;
      this->topPriorityAnnotationTypeFreq.at(i, &t, &freq);
      fprintf(fp, "%s\t%d\n", AnnotationString[t], freq);
    }
    fclose(fp);
  };
  
  void printBaseChangeFrequency(const char* fileName){
    FILE* fp = fopen(fileName, "wt");
    assert(fp);
    unsigned int n = this->baseFreq.size();
    for (unsigned int i = 0; i < n; i++){
      std::string k;
      int freq;
      this->baseFreq.at(i, &k, &freq);
      fprintf(fp, "%s\t%d\n", k.c_str(), freq);
    }
    fclose(fp);
  };
  void printCodonChangeFrequency(const char* fileName){
    FILE* fp = fopen(fileName, "wt");
    assert(fp);
    unsigned int n = this->codonFreq.size();
    for (unsigned int i = 0; i < n; i++){
      std::string k;
      int freq;
      this->codonFreq.at(i, &k, &freq);
      fprintf(fp, "%s\t%d\n", k.c_str(), freq);
    }
    fclose(fp);
  };
  void printIndelLengthFrequency(const char* fileName){
    FILE* fp = fopen(fileName, "wt");
    assert(fp);
    unsigned int n = this->indelLengthFreq.size();
    for (unsigned int i = 0; i < n; i++){
      int key;
      int freq;
      this->indelLengthFreq.at(i, &key, &freq);
      fprintf(fp, "%d\t%d\n", key, freq);
    }
    fclose(fp);
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
  void findInRangeGene(const std::string& chr, const int pos, std::vector<unsigned int>* potentialGeneIdx) {
    assert(potentialGeneIdx);
    potentialGeneIdx->clear();

    std::vector<Gene>& g = this->geneList[chr];
    unsigned int gLen = g.size();
    if (gLen == 0) {
      return;
    }
    int maxDist = (param.upstreamRange > param.downstreamRange) ? param.upstreamRange : param.downstreamRange;
    Range r ((pos - maxDist), (pos + maxDist));
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
    const Chromosome& seq = this->gs[chr];
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
  void annotateIns(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt, AnnotationResult* result) {
    Gene& g = this->geneList[chr][geneIdx];
    result->add(g);

    // might useful vars.
    int dist2Gene;
    int utrPos, utrLen;
    int exonNum; // which exon
    int codonNum; // which codon
    int codonPos[3] = {0, 0, 0}; // the codon position
    AnnotationType type; // could be one of
    int intronNum; // which intron
    bool isEssentialSpliceSite;

    result->add(INSERTION);
    if (g.isUpstream(variantPos, param.upstreamRange, &dist2Gene)) {
      result->add(UPSTREAM);
    } else if (g.isDownstream(variantPos, param.upstreamRange, &dist2Gene)) {
      result->add(DOWNSTREAM);
    } else if (g.isExon(variantPos, &exonNum)){//, &codonNum, codonPos)) {
      result->add(EXON);
      if (g.isCoding()) {
        if (g.is5PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR5);
        } else if (g.is3PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR3);
        } else { // cds part has base change
          int insertSize = alt.size() - ref.size();
          if (insertSize % 3 == 0) {
            result->add(CODON_GAIN);
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
            result->addDetail(CODON_GAIN, s);
          } else {
            result->add(FRAME_SHIFT);
          }
        }
      } else {
        result->add(NONCODING);
      }
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
      }
    } else if (g.isIntron(variantPos, &intronNum)) {
      result->add(INTRON);
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
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
  void annotateDel(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt, AnnotationResult* result) {
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
    result->add(g);
    for (std::set<AnnotationType>::const_iterator it = annotationSet.begin();
         it != annotationSet.end();
         it++) {
      result->add(*it);
    };
  }; // end annotateDel
  /**
   * SV is the most complex scenario. fully support this is an ongoing work.
   * We will just annotation the region in the rough scale.
   */
  void annotateSV(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt, AnnotationResult* result) {
    Gene& g = this->geneList[chr][geneIdx];
    result->add(g);

    // might useful vars.
    int dist2Gene;
    int utrPos, utrLen;
    int exonNum; // which exon
    int codonNum; // which codon
    int codonPos[3] = {0, 0, 0}; // the codon position
    AnnotationType type; // could be one of
    int intronNum; // which intron
    bool isEssentialSpliceSite;
    result->add(STRUCTURE_VARIATION);
    if (g.isUpstream(variantPos, param.upstreamRange, &dist2Gene)) {
      result->add(UPSTREAM);
    } else if (g.isDownstream(variantPos, param.upstreamRange, &dist2Gene)) {
      result->add(DOWNSTREAM);
    } else if (g.isExon(variantPos, &exonNum)){//, &codonNum, codonPos)) {
      result->add(EXON);
      if (!g.isNonCoding()) {
        if (g.is5PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR5);
        } else if (g.is3PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR3);
        } else { // cds part has base change
          result->add(CODON_REGION);
        }
      } else{
        result->add(NONCODING);
      }
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
      }
    } else if (g.isIntron(variantPos, &intronNum)) {
      result->add(INTRON);
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
      }
    } else {
      //annotation.push_back("Intergenic");
    }
  } // end annotateSV
  void annotateSNP(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt, AnnotationResult* result) {
    Gene& g = this->geneList[chr][geneIdx];
    result->add(g);

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
      result->add(UPSTREAM);
    } else if (g.isDownstream(variantPos, param.upstreamRange, &dist2Gene)) {
      result->add(DOWNSTREAM);
    } else if (g.isExon(variantPos, &exonNum)){//, &codonNum, codonPos)) {
      result->add(EXON);
      if (g.isCoding()) {
        if (g.is5PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR5);
        } else if (g.is3PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR3);
        } else { // cds part has base change
          if (g.calculateCodonPosition(variantPos, &codonNum, codonPos)) {
            char refTriplet[3];
            char altTriplet[3];
            std::string refAAName;
            std::string altAAName;
            std::string refLetterName;
            std::string altLetterName;
            AnnotationType annotationType;
            // when reference genome is provided
            if (this->gs.size() >0){
              this->fillTriplet(chr, variantPos, codonPos, g.forwardStrand, ref, alt, refTriplet, altTriplet);
              if (!g.forwardStrand){
                complementTriplet(refTriplet);
                complementTriplet(altTriplet);
              };
              refAAName = this->codon.toAA(refTriplet);
              altAAName = this->codon.toAA(altTriplet);
              refLetterName = this->codon.toLetter(refTriplet);
              altLetterName = this->codon.toLetter(altTriplet);
              annotationType = this->determineSNVType(refAAName, altAAName, codonNum);

              result->add(annotationType);
              std::string s;
              s += WITHIN_GENE_LEFT_DELIM;
              s += refTriplet[0];
              s += refTriplet[1];
              s += refTriplet[2];
              s += "/";
              s += refAAName;
              s += "/";
              s += refLetterName;
              s += "->";
              s += altTriplet[0];
              s += altTriplet[1];
              s += altTriplet[2];
              s += "/";
              s += altAAName;
              s += "/";
              s += altLetterName;
              s += WITHIN_GENE_SEPARATOR;
              // quick patch about codon number
              char buf[128];
              sprintf(buf, "Base%d/%d",
                      codonNum + 1, g.getCDSLength());
              s += buf;
              s += WITHIN_GENE_SEPARATOR;
              s += "Codon";
              s += toStr( codonNum / 3 + 1 );
              s += "/";
              s += toStr(g.getCDSLength() / 3);
              s += WITHIN_GENE_SEPARATOR;
              s += "Exon";
              s += toStr(exonNum + 1); // convert 0 indexed to 1 indexed
              s += "/";
              s += toStr( (int)( g.exon.size()));
              s += WITHIN_GENE_RIGHT_DELIM;
              result->addDetail(annotationType, s);
              // record frequency
              this->codonFreq.add(refAAName+"->"+altAAName);
            } else {
              result->add(SNV);
            }
          } else {
            result->add(SNV);
          }
        }
      } else{
        result->add(NONCODING);
      }
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
      }
    } else if (g.isIntron(variantPos, &intronNum)) {
      result->add(INTRON);
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
      }
    } else {
      //annotation.push_back("Intergenic");
    }
  } // end annotateSNP
  /**
   * annotation results will be store in @param result
   */
  void annotateByGene(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt,
                      const VARIATION_TYPE& type,
                      AnnotationResult* result){
    result->clear();
    switch(type) {
      case SNP:
        this->annotateSNP(geneIdx, chr, variantPos, ref, alt, result);
        break;
      case INS:
        this->annotateIns(geneIdx, chr, variantPos, ref, alt, result);
        break;
      case DEL:
        this->annotateDel(geneIdx, chr, variantPos, ref, alt, result);
        break;
      case SV:
        this->annotateSV(geneIdx, chr, variantPos, ref, alt, result);
        break;
      case MIXED:
      case UNKNOWN:
      default:
        LOG << "Currently we don't support this variation type: " << type << "\n";
        break;
    };

    return;
  };
  /**
   * @return indel length, for insertion, return positive number; or return negative number
   */
  int calculateIndelLength(const std::string& ref, const std::string& alt){
    int refLen = ref.size();
    int altLen = alt.size();
    if (alt == "." || alt == "<DEL>") {
      altLen = 0 ;
    }
    return (altLen - refLen);
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
  Priority priority;
  GeneFormat format;
  bool allowMixedVariation;       // VCF ALT field may have more than one variation e..g A,C

  std::vector<AnnotationResult> annotationResult;
  AnnotationOutput outputter;               // control output format

  FreqTable<AnnotationType> annotationTypeFreq;                 // base change frequency
  FreqTable<AnnotationType> topPriorityAnnotationTypeFreq;      // base change frequency of top priority
  FreqTable<std::string> baseFreq;                              // base change frequency
  FreqTable<std::string> codonFreq;                             // codon change frequency
  FreqTable<int> indelLengthFreq;                               // for insertion, the value is positive; for deletion, positive

}; // end class GeneAnnotation

int main(int argc, char *argv[])
{
  banner(stdout);

  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Required Parameters")
      ADD_STRING_PARAMETER(pl, inputFile, "-i", "Specify input VCF file")
      ADD_STRING_PARAMETER(pl, outputFile, "-o", "Specify output VCF file")
      ADD_STRING_PARAMETER(pl, geneFile, "-g", "Specify gene file")
      ADD_PARAMETER_GROUP(pl, "Optional Parameters")
      ADD_STRING_PARAMETER(pl, inputFormat, "--inputFormat", "Specify format (default: vcf). \"-f plain \" will use first 4 columns as chrom, pos, ref, alt")
      ADD_STRING_PARAMETER(pl, referenceFile, "-r", "Specify reference genome position")
      ADD_STRING_PARAMETER(pl, geneFileFormat, "-f", "Specify gene file format (default: refFlat, other options knownGene)")
      ADD_STRING_PARAMETER(pl, priorityFile, "-p", "Specify priority of annotations")
      ADD_STRING_PARAMETER(pl, codonFile, "-c", "Specify codon file (default: codon.txt)")
      ADD_INT_PARAMETER(pl, upstreamRange, "-u", "Specify upstream range (default: 50)")
      ADD_INT_PARAMETER(pl, downstreamRange, "-d", "Specify downstream range (default: 50)")
      ADD_INT_PARAMETER(pl, spliceIntoExon, "--se", "Specify splice into extron range (default: 3)")
      ADD_INT_PARAMETER(pl, spliceIntoIntron, "--si", "Specify splice into intron range (default: 8)")
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
  GeneAnnotationParam param;
  param.upstreamRange = FLAG_upstreamRange ? FLAG_upstreamRange : 50;
  param.downstreamRange = FLAG_downstreamRange ? FLAG_downstreamRange : 50;
  param.spliceIntoExon = FLAG_spliceIntoExon ? FLAG_spliceIntoExon : 3;
  param.spliceIntoIntron = FLAG_spliceIntoIntron ? FLAG_spliceIntoIntron : 8;

  std::string logFileName = FLAG_outputFile + ".log";
  LOG_START(logFileName.c_str());
  LOG_START_TIME;
  LOG_PARAMETER(pl);

  GeneAnnotation ga;
  pl.Status();
  ga.setAnnotationParameter(param);

  if (FLAG_priorityFile.size() == 0) {
    fprintf(stderr, "Use default priority file: /net/fantasia/home/zhanxw/anno/priority.txt\n");
    FLAG_priorityFile = "/net/fantasia/home/zhanxw/anno/priority.txt";
  };
  if (FLAG_codonFile.size() == 0) {
    fprintf(stderr, "Use default codon file: /net/fantasia/home/zhanxw/anno/codon.txt\n");
    FLAG_codonFile = "/net/fantasia/home/zhanxw/anno/codon.txt";
  }
  if (FLAG_referenceFile.size() == 0) {
    fprintf(stderr, "Use default codon file: /net/fantasia/home/zhanxw/anno/codon.txt\n");
    FLAG_referenceFile = "/data/local/ref/karma.ref/human.g1k.v37.fa";
  }

  ga.openReferenceGenome(FLAG_referenceFile.c_str());
  ga.openCodonFile(FLAG_codonFile.c_str());
  ga.openPriorityFile(FLAG_priorityFile.c_str());

  ga.setFormat(FLAG_geneFileFormat);
  ga.openGeneFile(FLAG_geneFile.c_str());
  if (toLower(FLAG_inputFormat) == "VCF" || FLAG_inputFormat.size() == 0) {
    ga.annotateVCF(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
  } else if (toLower(FLAG_inputFormat) == "PLAIN") {
    ga.annotatePlain(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
  } else if (toLower(FLAG_inputFormat) == "PLINK") {
    ga.annotatePlink(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
  } else{
    fprintf(stderr, "Cannot recognize input file format: %s \n", FLAG_inputFile.c_str());
    abort();
  };

  LOG_END_TIME;
  LOG_END ;
  printf("Annotation succeed!\n");
  return 0;
}
