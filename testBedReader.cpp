#include "IO.h"
#include "TypeConversion.h"
#include "StringUtil.h"
#include <map>
#include <string>
#include <vector>
#include <algorithm>

#include "BedReader.h"

int main(int argc, char *argv[])
{
  BedReader br;
  int r = br.open("test.bed");
  printf("load %d regions\n", r);
  // br.dump();
  std::vector<std::string> ret;
  if (br.find("1", 100, &ret)) {
    printf("found 1:100 at %s\n", stringJoin(ret, ',').c_str());
  } else {
    printf("not found 1:100\n");
  }
  
  return 0;
}
