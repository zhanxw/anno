#include "GenomeScore.h"
#include "TypeConversion.h"

int main(int argc, char *argv[])
{
  const char path[] = "/net/fantasia/home/zhanxw/anno/resources/gerp/";
  printf("Open genome score file: %s\n", path);
  GenomeScore gs("/net/fantasia/home/zhanxw/anno/resources/gerp/");
  //float baseScore(const char* chrom, int pos) {
  for (int i = 1; i <= 10; ++i){
    printf("1 %d %f\n", i, gs.baseScore("1", i));
  }
  const int chr1Len = 249250621;
  for (int i = 0; i <= 10; ++i){
    printf("1 %d %f\n", chr1Len - i, gs.baseScore("1", chr1Len - i));
  }
  
  for (int i = 1; i <= 10; ++i){
    printf("22 %d %f\n", i, gs.baseScore("22", i));
  }
  const int chr22Len = 51304566;
  for (int i = 0; i <= 10; ++i){
    printf("22 %d %f\n", chr22Len - i , gs.baseScore("22", chr22Len - i));
  }

  int pos;
  pos = 110199012;
  printf("1 %d %f\n", pos, gs.baseScore("1", pos));
  pos = 110199013;
  printf("1 %d %f\n", pos, gs.baseScore("1", pos));
  
  return 0;
}
