#ifndef _ANNOTATIONINPUTFILE_H_
#define _ANNOTATIONINPUTFILE_H_


typedef enum {VCF = 0 , PLAIN, PLINK, EPACTS} InputFileFormat;

// hold input file
// and return (chrom, pos, ref, alt) iteratively
class AnnotationInputFile{
 public:
  AnnotationInputFile(const char* inputFileName, const char* inputFormatStr) {
    // check inputFormat
    std::string inputFormat = toLower(inputFormatStr);
    if (!inputFormat.empty() &&
        inputFormat != "vcf" &&
        inputFormat != "plain" &&
        inputFormat != "plink" &&
        inputFormat != "epacts") {
      fprintf(stderr, "Unsupported input format [ %s ], we support VCF, plain, plink and EPACTS formats.\n", inputFormatStr);
      LOG << "Unsupported input format [ " << inputFormatStr << " ], we support VCF, plain, plink and EPACTS formats.\n";
      abort();
    };

    if (inputFormat == "vcf" || inputFormat.empty()) {
      // ga.annotateVCF(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = VCF;
    } else if (inputFormat == "plain") {
      // ga.annotatePlain(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = PLAIN;
    } else if (inputFormat == "plink") {
      // ga.annotatePlink(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = PLINK;
    } else if (inputFormat == "epacts") {
      // ga.annotateEpacts(FLAG_inputFile.c_str(), FLAG_outputFile.c_str());
      this->format = EPACTS;
    } else{
      fprintf(stderr, "Cannot recognize input file format: %s \n", inputFileName);
      abort();
    };

    // open input files
    this->lr = new LineReader(inputFileName);

    // init values
    checkReference = true;
    failedReferenceSite = 0;
  };
  ~AnnotationInputFile() {
    this->close();
  }
  void close() {
    if (this->lr) {
      delete lr;
      this->lr = NULL;
    }
    if (checkReference && failedReferenceSite > 0) {
      fprintf(stderr, "ERROR: Total [ %d ] sites have unmatched Reference alleles\n", failedReferenceSite);
      LOG << "ERROR: Total [ " << failedReferenceSite << " ] sites have unmatched Reference alleles\n";
    }
  };


  int openReferenceGenome(const char* referenceGenomeFileName) {
    return this->gs.open(referenceGenomeFileName);
  }
  // // extract one header line
  // bool extractHeader(std::string* l) {
  //   this->line.clear();
  //   bool ret = this->lr->readByLine(this->line);
  //   *l = this->line;
  //   return ret;
  // };
  // // extract one header line
  // bool extractHeader(std::string* l) {
  //   this->line.clear();
  //   bool ret = this->lr->readByLine(this->line);
  //   *l = this->line;
  //   return ret;
  // };

  // check ref and alt alleles (may switch ref and alt)
  //@return true if (1) ref match reference; (2) after switch ref and alt, match reference genome.
  bool forceReferenceStrand(const std::string& chrom,
                            const int& pos,
                            std::string* ref,
                            std::string* alt) {
    // determine ref base from reference
    bool refMatchRef = true;
    for (size_t i = 0; i < ref->size(); i++ ) {
      if ((*ref)[i] != gs[chrom][pos - 1 + i]) {
        refMatchRef = false;
        break;
      }
    }
    if (!refMatchRef) {
      bool altMatchRef = true;
      for (size_t i = 0; i < alt->size(); i++ ) {
        if ( (*alt)[i] != gs[chrom][pos - 1 + i]) {
          altMatchRef = false;
          break;
        }
      }
      if (!altMatchRef) {
        fprintf(stderr, "Ref [%s] and alt [%s] does not match reference: %s:%d\n", ref->c_str(), alt->c_str(),
                chrom.c_str(), pos);
        return false;
      } else {
        std::swap(*ref, *alt);
      }
    }
    return true;
  }
  void setCheckReference(bool b) {
    this->checkReference = b;
  };
  // if reach end or experience something wrong, will @return false
  bool extract(std::string* chrom,
               int* pos,
               std::string* ref,
               std::string* alt) {
    bool ret;
    do {
      ret = this->lr->readLine(&this->line);
      if (ret == false)
        return ret;
    } while (this->line.empty());

    // for any line beginning with '#', store headers
    // " CHR " is for PLINK header, "CHROM" is for other headers    
    while (line[0] == '#' || line.substr(0,5) == "CHROM" || line.substr(0, 5) == " CHR ") { 
      this->header.push_back(line);
      do {
        ret = this->lr->readLine(&this->line);
        if (ret == false)
          return ret;
      } while (this->line.empty());
    }



    switch (this->format){
      case VCF:
        stringTokenize(line, "\t ", &fd);
        if (fd.size() < 5) return false;
        *chrom = fd[0];
        *pos = toInt(fd[1]);
        *ref = fd[3];
        *alt = fd[4];
        break;
      case PLAIN:
        stringNaturalTokenize(line, "\t ", &fd);
        if (fd.size() < 4) return false;
        *chrom = fd[0];
        *pos = toInt(fd[1]);
        *ref = fd[2];
        *alt = fd[3];
        break;
      case PLINK:
        stringNaturalTokenize(line, "\t ", &fd);
        if (fd.size() < 10) return false;
        *chrom = fd[0];
        *pos = toInt(fd[2]);
        *ref = fd[3];
        *alt = fd[6];

        if (!forceReferenceStrand(*chrom, *pos, ref, alt))
          return false;
        break;
      case EPACTS:
        {
          stringNaturalTokenize(line, "\t ", &fd);
          // e.g.
          // 20:139681_G/A   266     1       1       0.0018797       NA      NA
          // find
          int beg = 0;
          int sep = fd[0].find(':', beg);
          *chrom = fd[0].substr(beg, sep - beg);

          beg = sep + 1;
          sep = fd[0].find('_', beg);
          *pos = toInt(fd[0].substr(beg, sep - beg));

          beg = sep + 1;
          sep = fd[0].find('/', beg);
          *ref = toUpper(fd[0].substr(beg, sep - beg));

          beg = sep + 1;
          sep = fd[0].find_first_of(" _\t", beg);
          *alt = toUpper(fd[0].substr(beg, sep - beg));

          epactsPrefixLength = sep;
          if ( chrom->empty() || *pos <= 0 || ref->empty() || alt->empty()) {
            fprintf(stderr, "Skip line: %s ...." , fd[0].c_str());
            LOG << "Skip: " << fd[0];
            return false;
          }
        }
        break;
      default:
        fprintf(stderr, "Unknown format, quitting!\n");
        abort();
        break;
    }// end switch

    // verify reference
    if (this->checkReference) {
      std::string refFromGenome = this->gs.getBase(*chrom, *pos, *pos + ref->size());
      if ((*ref) != refFromGenome) {
        ++ failedReferenceSite;
        if (failedReferenceSite <= 10) { // output at most 10 warnings
          fprintf(stderr, "ERROR: Reference allele does not match genome reference [ %s:%d %s]\n", chrom->c_str(), *pos, ref->c_str());
        }
        LOG << "ERRROR: Reference allele [" << ref <<   "]  does not match reference genome [" << refFromGenome << "] at " << *chrom << ":" << *pos << "\n";
      };
    }
    return true;
  };
  InputFileFormat getFormat() const {
    return this->format;
  };
  std::string getEpactsPrefix(const std::string& s) const{
    return s.substr(0, this->epactsPrefixLength);
  };
  const std::vector< std::string>& getHeader() const {
    return this->header;
  };
  const std::vector< std::string>& getFields() const {
    return this->fd;
  };
 private:
  bool checkReference;
  int failedReferenceSite;
  InputFileFormat format;
  LineReader* lr;
  std::vector< std::string> fd;
  std::string line;
  std::vector< std::string> header;
  GenomeSequence gs; // check if ref alleles matched reference genome
  size_t epactsPrefixLength;
}; // end class AnnotationInputFile



#endif /* _ANNOTATIONINPUTFILE_H_ */
