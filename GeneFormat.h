#ifndef _GENEFORMAT_H_
#define _GENEFORMAT_H_
//TODO
//Use format.txt to flexible specify inputs
class GeneFormat{
public:
  void setRefFlatFormat(){
    this->nameCol.clear();
    this->nameCol.push_back(0);
    if (!this->skipTranscriptName)
      this->nameCol.push_back(1);
    this->chrCol = 2;
    this->strandCol = 3;
    this->txStartCol = 4;
    this->txEndCol = 5;
    this->cdsStartCol = 6;
    this->cdsEndCol = 7;
    this->exonNumCol = 8;
    this->exonStartCol = 9;
    this->exonEndCol = 10;

    this->expectedColumn = 11;
  };
  void setUCSCKnownGeneFormat(){
    this->nameCol.clear();
    this->nameCol.push_back(10);
    if (!this->skipTranscriptName)
      this->nameCol.push_back(1);
    this->nameCol.push_back(0);
    this->chrCol = 1;
    this->strandCol = 2;
    this->txStartCol = 3;
    this->txEndCol = 4;
    this->cdsStartCol = 5;
    this->cdsEndCol = 6;
    this->exonNumCol = 7;
    this->exonStartCol = 8;
    this->exonEndCol = 9;

    this->expectedColumn = 11;
  };
  void setRefGeneFormat() {
    this->nameCol.clear();
    this->nameCol.push_back(12);
    if (!this->skipTranscriptName)
      this->nameCol.push_back(1);
    this->nameCol.push_back(1);
    this->chrCol = 2;
    this->strandCol = 3;
    this->txStartCol = 4;
    this->txEndCol = 5;
    this->cdsStartCol = 6;
    this->cdsEndCol = 7;
    this->exonNumCol = 8;
    this->exonStartCol = 9;
    this->exonEndCol = 10;

    this->expectedColumn = 16;
  }
  int getExpectedColumn() const {
    return this->expectedColumn;
  };
  /**
   * Optionally hide or show transcript names
   * NOTE: call these before set any other formats.
   */
  void hideTranscriptName() {
    this->skipTranscriptName = true;
  };
  void showTranscriptName() {
    this->skipTranscriptName = false;
  };
  
public:
  std::vector<int> nameCol;
  int chrCol;
  int strandCol;
  int txStartCol;
  int txEndCol;
  int cdsStartCol;
  int cdsEndCol;
  int exonNumCol;
  int exonStartCol;
  int exonEndCol;
  int expectedColumn;

  bool skipTranscriptName;
};

#endif /* _GENEFORMAT_H_ */
