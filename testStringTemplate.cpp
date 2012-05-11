#include "StringTemplate.h"
#include <cassert>

int main(int argc, char *argv[])
{
  std::string ret;
  StringTemplate t1("");
  std::string e1 = "";
  assert( 0 == t1.translate(&ret));
  assert( ret == e1);

  t1 = "a";
  e1 = "a";
  assert( 0 == t1.translate(&ret));
  assert( ret == e1);

  t1 = "01$(a)34";
  e1 = "01234";
  t1.add("a", "2");
  assert( 0 == t1.translate(&ret));
  assert( ret == e1);


  t1 = "01$[a]34";
  e1 = "01a34";
  assert( 0 == t1.translate(&ret));
  assert( ret == e1);


  t1 = "01$[23$(a)56]78";
  e1 = "012345678";
  t1.add("a", "4");
  assert( 0 == t1.translate(&ret));
  assert( ret == e1);


  t1 = "01$[$(a)]45";
  e1 = "012345";
  std::vector<std::string> a;
  a.push_back("2");
  a.push_back("3");
  t1.add("a", a);
  assert( 0 == t1.translate(&ret));
  assert( ret == e1);


  t1 = "0,1,$[$(a) ,],4,5";
  e1 = "0,1,2,3,4,5";
  a.clear();
  a.push_back("2");
  a.push_back("3");
  t1.add("a", a);
  assert( 0 == t1.translate(&ret));
  assert( ret == e1);

  t1 = "$[$(a) ,],$[$(b) ,],$(c)";
  e1 = "0,1,2,3,4,5";
  a.clear();
  a.push_back("0");
  a.push_back("1");
  t1.add("a", a);
  std::vector<std::string> b;
  b.push_back("2");
  b.push_back("3");
  b.push_back("4");
  t1.add("b", b);
  t1.add("c", "5");
  assert( 0 == t1.translate(&ret));
  assert( ret == e1);

  
  return 0;
}
