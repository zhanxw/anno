#include <stdio.h>
#include "LineBreaker.h"

int main(int argc, char *argv[])
{
  std::string s = "aaa bb cc ddddd";
  LineBreaker lb(6);
  lb.setContent(s);
  size_t h = lb.getHeight();
  for (size_t i = 0; i < h; ++i) {
    printf("%s\n", lb[i].c_str());
  }
  return 0;
}
