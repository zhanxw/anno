/**
 * 3. Output in VAT format
 */
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

#include "Type.h"
#include "Codon.h"
#include "GenomeSequence.h"
#include "Gene.h"
#include "GeneFormat.h"
#include "SequenceUtil.h"
#include "FreqTable.h"
#include "StringTemplate.h"
#include "GeneAnnotation.h"
#include "BedReader.h"
#include "GenomeScore.h"
#include "ModelParser.h"
#include "TabixReader.h"
#include "GitVersion.h"

#include "AnnotationOutputFile.h"

void banner(FILE* fp) {
  const char* string =
      "..............................................         \n"
      " ...      Anno(tation)                       ...       \n"
      "  ...      Xiaowei Zhan, Goncalo Abecasis     ...      \n"
      "   ...      Speical Thanks:                    ...     \n"
      "    ...      Hyun Ming Kang, Yanming Li         ...    \n"
      "     ...      zhanxw@umich.edu                    ...  \n"
      "      ...      Sep 2011                            ... \n"
      "       ................................................\n"
      "                                                       \n"
      ;
  fputs(string, fp);
#ifndef NDEBUG
  const char* debug =
      "-------------------------------------------------------\n"
      "|                                                     |\n"
      "|                DEBUG  MODE                          |\n"
      "|                (slow mode)                          |\n"
      "|                                                     |\n"
      "|   Try:                                              |\n"
      "|      make clean; make release                       |\n"
      "|   then run:                                         |\n"
      "|      ./executable/anno                              |\n"
      "|   to use release/fast version                       |\n"
      "|                                                     |\n"
      "-------------------------------------------------------\n"
      ;
  fputs(debug, fp);
#endif
};

extern const char* gitVersion;

OutputAnnotationString AnnotationString; // global variable


// run annotations (gene based, bed file, genomeScores, tabix database)
// and store results
class AnnotationController{
 public:
  AnnotationController(AnnotationInputFile& in): aif(in){
  };
  virtual ~AnnotationController() {
    for (size_t i = 0; i < bedReader.size() ; ++i) {
      delete bedReader[i];
    }
    for (size_t i = 0; i < genomeScore.size(); ++i ) {
      delete genomeScore[i];
    };
    for (size_t i = 0; i < tabixReader.size(); ++i ) {
      delete tabixReader[i];
    };

    
  };
  void openBedFile(const char* tag, const char* fn) {
    // check duplication
    for (size_t i = 0; i < this->bedTag.size(); ++i ) {
      if (this->bedTag[i] == tag) {
        fprintf(stderr, "ERROR: Duplicated tag [ %s ] \n", tag);
        return;
      }
    }

    // add bedFile
    BedReader* p = new BedReader;
    int ret = p->open(fn);
    if (ret < 0) {
      fprintf(stderr, "Cannot open BED file: [ %s ]\n", fn);
      delete p;
      return ;
    } else {
      fprintf(stderr, "DONE: Load %d regions from BED file\n", ret);
    };

    this->bedTag.push_back(tag);
    this->bedReader.push_back(p);
  };
  void openGenomeScoreFile(const char* tag, const char* fn){
    // check duplication
    for (size_t i = 0; i < this->genomeScoreTag.size(); ++i) {
      if (this->genomeScoreTag[i] == tag ) {
        fprintf(stderr, "ERROR: Duplicated tag [ %s ] \n", tag);
        return;
      }
    }

    // add genome score file
    GenomeScore* p = new GenomeScore(fn);
    this->genomeScoreTag.push_back(tag);
    this->genomeScore.push_back(p);
  };

  void addTabixReader(TabixReader* t) {
    this->tabixReader.push_back(t);
  };

  void annotate(std::string& chrom,
                int& pos,
                std::string& ref,
                std::string& alt) {
    this->geneAnnotation.annotate(chrom, pos, ref, alt);

    this->result.clear();
    this->result["ANNO"] = this->geneAnnotation.getTopPriorityAnnotation();
    this->result["ANNOFULL"] = this->geneAnnotation.getFullAnnotation();

    std::vector<std::string> bedString;
    for (size_t br = 0; br < this->bedReader.size(); ++ br) {
      if (this->bedReader[br]->find(chrom.c_str(), pos, &bedString)){
        if (!bedString.empty())  {
          this->result[this->bedTag[br]] = stringJoin(bedString, ',');
        } else {
          this->result[this->bedTag[br]] = "";
        }
      }
    }
    for (size_t gs = 0; gs < this->genomeScore.size(); ++ gs) {
      this->result[this->genomeScoreTag[gs]] = toStr(this->genomeScore[gs]->baseScore(chrom.c_str(), pos));
    }
    for (size_t tb = 0; tb < this->tabixReader.size(); ++ tb) {
      TabixReader& tabix = *this->tabixReader[tb];
      tabix.addAnnotation(chrom, pos, ref, alt);
      size_t s = tabix.getTag().size();
      for (size_t i = 0; i < s; ++i) {
        this->result[tabix.getTag()[i]] = tabix.getAnnotation()[i];
      }
    }
  };

  const OrderedMap<std::string, std::string>& getResult() const {
    return this->result;
  }
 public:
  AnnotationInputFile& aif;
  // various annotation types
  GeneAnnotation geneAnnotation;
 private:
  // various annotation types
  std::vector<BedReader*> bedReader;
  std::vector<std::string> bedTag;
  std::vector<GenomeScore*> genomeScore;
  std::vector<std::string> genomeScoreTag;
  std::vector<TabixReader*> tabixReader;
  OrderedMap<std::string, std::string> result; // store all types of annotation results
};

int main(int argc, char *argv[])
{
  banner(stderr);

  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Required Parameters")
      ADD_STRING_PARAMETER(pl, inputFile, "-i", "Specify input file")
      ADD_STRING_PARAMETER(pl, outputFile, "-o", "Specify output file")
      ADD_PARAMETER_GROUP(pl, "Gene Annotations")
      ADD_STRING_PARAMETER(pl, geneFile, "-g", "Specify gene file")
      ADD_STRING_PARAMETER(pl, referenceFile, "-r", "Specify reference genome position")      
      ADD_STRING_PARAMETER(pl, inputFormat, "--inputFormat", "Specify format (default: vcf). \"-f plain \" will use first 4 columns as chrom, pos, ref, alt")
      ADD_BOOL_PARAMETER(pl, checkReference, "--checkReference", "Check whether reference alleles matches genome reference")
      ADD_STRING_PARAMETER(pl, geneFileFormat, "-f", "Specify gene file format (default: refFlat, other options: knownGene, refGene)")
      ADD_STRING_PARAMETER(pl, priorityFile, "-p", "Specify priority of annotations")
      ADD_STRING_PARAMETER(pl, codonFile, "-c", "Specify codon file (default: codon.txt)")
      ADD_INT_PARAMETER(pl, upstreamRange, "-u", "Specify upstream range (default: 50)")
      ADD_INT_PARAMETER(pl, downstreamRange, "-d", "Specify downstream range (default: 50)")
      ADD_INT_PARAMETER(pl, spliceIntoExon, "--se", "Specify splice into extron range (default: 3)")
      ADD_INT_PARAMETER(pl, spliceIntoIntron, "--si", "Specify splice into intron range (default: 8)")
      ADD_PARAMETER_GROUP(pl, "Other Annotations")
      ADD_STRING_PARAMETER(pl, genomeScore, "--genomeScore", "Specify the folder of genome score (e.g. GERP=dirGerp/,SIFT=dirSift/)")
      ADD_STRING_PARAMETER(pl, bedFile, "--bed", "Specify the bed file and tag (e.g. ONTARGET1=a1.bed,ONTARGET2=a2.bed)")
      ADD_STRING_PARAMETER(pl, tabixFile, "--tabix", "Specify the tabix file and tag (e.g. abc.txt.gz(chrom=1,pos=7,ref=3,alt=4,SIFT=7,PolyPhen=10)")      
      ADD_PARAMETER_GROUP(pl, "Auxillary Functions")      
      ADD_STRING_PARAMETER(pl, outputFormat, "--outputFormat", "Specify predefined annotation words (default or epact)")
      ADD_BOOL_PARAMETER(pl, indexOutput, "--indexOutput", "Specify whether to index output file using tabix (require .gz suffix for output file)")      
      END_PARAMETER_LIST(pl)
      ;

  pl.Read(argc, argv);
  if (FLAG_geneFileFormat.empty()) {
    FLAG_geneFileFormat = "refFlat";
  }

  if (FLAG_inputFile.empty()) {
    pl.Help();
    fprintf(stderr, "Please specify input file\n");
    exit(1);
  }
  if (FLAG_outputFile.empty()) {
    pl.Help();
    fprintf(stderr, "Please specify output file\n");
    exit(1);
  }
  if (!hasSuffix(FLAG_outputFile,".gz") && FLAG_indexOutput) {
    fprintf(stderr, "Please give output file \".gz\" suffix to enable index (--indexOutput).\n");
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
  LOG << "Version: " << gitVersion << "\n";
  
  pl.Status();
  if (FLAG_REMAIN_ARG.size() > 0){
    fprintf(stderr, "Unparsed arguments: ");
    for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++){
      fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
    }
    abort();
  }
  
  if (FLAG_geneFile.empty()) {
    pl.Help();
    fprintf(stderr, "Please specify gene file\n");
    exit(1);
  }
  
  if (FLAG_priorityFile.empty()) {
    fprintf(stderr, "Use default priority file: /net/fantasia/home/zhanxw/anno/priority.txt\n");
    FLAG_priorityFile = "/net/fantasia/home/zhanxw/anno/priority.txt";
  };
  if (FLAG_codonFile.empty()) {
    fprintf(stderr, "Use default codon file: /net/fantasia/home/zhanxw/anno/codon.txt\n");
    FLAG_codonFile = "/net/fantasia/home/zhanxw/anno/codon.txt";
  }
  if (FLAG_referenceFile.empty()) {
    fprintf(stderr, "Use default codon file: /net/fantasia/home/zhanxw/anno/codon.txt\n");
    FLAG_referenceFile = "/data/local/ref/karma.ref/human.g1k.v37.fa";
  }

  if (!FLAG_outputFormat.empty()) {
    AnnotationString.setFormat(FLAG_outputFormat.c_str());
  }
  else {
    AnnotationString.setFormat("default");
  };


  AnnotationInputFile aif(FLAG_inputFile.c_str(), FLAG_inputFormat.c_str());
  aif.openReferenceGenome(FLAG_referenceFile.c_str());
  aif.setCheckReference(FLAG_checkReference);

  AnnotationController controller(aif);

  controller.geneAnnotation.setAnnotationParameter(param);
  controller.geneAnnotation.openReferenceGenome(FLAG_referenceFile.c_str());
  controller.geneAnnotation.openCodonFile(FLAG_codonFile.c_str());
  controller.geneAnnotation.openPriorityFile(FLAG_priorityFile.c_str());
  // controller.geneAnnotation.setFormat(FLAG_geneFileFormat);
  controller.geneAnnotation.openGeneFile(FLAG_geneFile.c_str(), FLAG_geneFileFormat.c_str());

  if (!FLAG_bedFile.empty()) {
    fprintf(stderr, "Use bed file: %s\n", FLAG_bedFile.c_str() );
    std::vector< std::string> fd;
    std::vector< std::string> bed;
    stringTokenize(FLAG_bedFile, ",", &fd);
    for (size_t i = 0; i < fd.size(); ++i ){
      stringTokenize(fd[i], "=", &bed);
      if (bed.size() == 2) {
        controller.openBedFile(bed[0].c_str(), bed[1].c_str());
      } else {
        fprintf(stderr, "ERROR: Cannot recognized format [ %s ].\n", fd[i].c_str());
        exit(1);
      };
    };
  };
  if (!FLAG_genomeScore.empty()){
    // fprintf(stderr, "Use binary GERP score: %s\n", FLAG_genomeScore.c_str());
    // ga.addGenomeScore("GERP", FLAG_genomeScore.c_str());
    fprintf(stderr, "Use binary score file: %s\n", FLAG_genomeScore.c_str() );
    std::vector< std::string> fd;
    std::vector< std::string> gs;
    stringTokenize(FLAG_genomeScore, ",", &fd);
    for (size_t i = 0; i < fd.size(); ++i ){
      stringTokenize(fd[i], "=", &gs);
      if (gs.size() == 2) {
        controller.openGenomeScoreFile(gs[0].c_str(), gs[1].c_str());
      } else {
        fprintf(stderr, "ERROR: Cannot recognized format [ %s ].\n", fd[i].c_str());
        exit(1);
      };
    };
  }
  // parse something like:
  // abc.txt.gz(chrom=1,pos=7,ref=3,alt=4,SIFT=7,PolyPhen=10)
  if(!FLAG_tabixFile.empty()){
    fprintf(stderr, "Use tabix file: %s\n", FLAG_tabixFile.c_str() );    
    ModelParser mp;
    mp.parse(FLAG_tabixFile);
    std::string fn = mp.getName();
    int chrom, pos, ref, alt;
    mp.assign("chrom", &chrom).assign("pos", &pos).assign("ref", &ref).assign("alt", &alt);
    fprintf(stderr, "Column %d, %d, %d and %d in tabix file will be matched to chromosome, position, reference allele, alternative allele respectively.\n", chrom, pos, ref, alt);
    TabixReader* tabix = new TabixReader(fn.c_str(), chrom, pos, ref, alt);
    
    for (size_t i = 0; i < mp.getParam().size(); ++i) {
      if ( toLower(mp.getParam()[i]) == "chrom" ||
           toLower(mp.getParam()[i]) == "pos" ||
           toLower(mp.getParam()[i]) == "ref" ||
           toLower(mp.getParam()[i]) == "alt") {
        continue;
      }
      int intValue;
      if (str2int(mp.getValue(i), &intValue)) {
        tabix->addTag(mp.getParam()[i], intValue);
        fprintf(stderr, "Tag %s will be from column %d in tabix file\n", mp.getParam()[i].c_str(), intValue);
      } else {
        tabix->addTag(mp.getParam()[i], mp.getValue(i));
        fprintf(stderr, "Tag %s will be from column %s (from header) in tabix file\n", mp.getParam()[i].c_str(), mp.getValue(i).c_str());        
      }
    }
    controller.addTabixReader(tabix);
  } // end halding tabix database

  std::string chrom;
  int pos;
  std::string ref;
  std::string alt;
  AnnotationOutputFile aof(FLAG_outputFile.c_str());
  aof.linkToInput(aif);
  std::string choppedChr; // chop leading 'chr'
  while (aif.extract(&chrom, &pos, &ref, &alt)) {
    choppedChr = chopChr(chrom);
    controller.annotate(choppedChr, pos, ref, alt);
    aof.writeResult(controller.getResult());
  };
  // aof.writeResult(controller.getResult()); // TODO: will add this to handle when input only have comment lines
  
  // output stats
  controller.geneAnnotation.outputAnnotationStats(FLAG_outputFile.c_str());

  aof.close();
  aif.close();
  LOG << "Annotate " << FLAG_inputFile << " to " << FLAG_outputFile << " succeed!\n";

  if (FLAG_indexOutput) {
    if (aof.indexOutput() == 0) {
      fprintf(stderr, "DONE: Indexing succeed!\n");
    } else {
      fprintf(stderr, "WARNING: Indexing failed!\n");
    }
  };

  LOG_END_TIME;
  LOG_END ;
  printf("Annotation succeed!\n");
  return 0;
}
