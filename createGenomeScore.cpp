#include "GenomeScore.h"
#include "TypeConversion.h"

void usage(int argc, char* argv[]) {
  printf("%s prefix suffix outputDir\n");
};

int main(int argc, char *argv[])
{
  std::string chrom;
  std::string prefix = argv[1]; //"/net/fantasia/home/zhanxw/software/gerp/chr";
  std::string suffix = argv[2]; //".maf.rates";
  std::string inFile;
  std::string outFile;
  for (int i = 1; i <= 22; ++i) {
    chrom = toString(i);
    inFile = prefix + chrom +suffix;
    outFile = argv[3];
    outFile += "chr";
    outFile += chrom;
    outFile += ".fbin";
    printf("Convert %s to %s\n", inFile.c_str(), outFile.c_str());
    GenomeScore::convert(chrom.c_str(), inFile.c_str(), outFile.c_str());
  }
  
  return 0;
}
