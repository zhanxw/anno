#include <cassert>
#include "TabixReader.h"

int main(int argc, char *argv[])
{
  {
    const char path[] = "/net/fantasia/home/zhanxw/anno/resources/dbNSFP/content/hg19/dbNSFP_light1.3.hg19.gz";    
    TabixReader tr(path, 1, 7, 3, 4); // dbNSFP chromosome

    tr.addAnnotation("1", 1, "A", "C");
    assert(tr.getAnnotation().size() == 0);

    tr.addTag("GerpNR", 16);
    tr.addTag("GerpRS", 17);
    tr.addTag("Last", 20);
    
    tr.addAnnotation("1", 1, "A", "C");
    assert(tr.getAnnotation()[0] == "NA");
    assert(tr.getAnnotation()[0] == tr.getAnnotation()[1]);
    assert(tr.getAnnotation()[0] == tr.getAnnotation()[2]);

    tr.addAnnotation("1", 69091, "A", "C");
    assert(tr.getAnnotation().size() == 3);
    printf("%s = %s\n", tr.getTag()[0].c_str(), tr.getAnnotation()[0].c_str());
    printf("%s = %s\n", tr.getTag()[1].c_str(), tr.getAnnotation()[1].c_str());
    printf("%s = %s\n", tr.getTag()[2].c_str(), tr.getAnnotation()[2].c_str());
  }

  {
    const char path[] = "/net/fantasia/home/zhanxw/anno/resources/dbNSFP/content/hg19/dbNSFP_light1.3.hg19.gz";    
    TabixReader tr(path, 1, 7, 3, 4); // dbNSFP chromosome
    tr.addTag("GerpNR", "GERP_NR");
    tr.addTag("GerpRS", "GERP_RS");

    tr.addAnnotation("1", 69091, "A", "C");
    assert(tr.getAnnotation().size() == 2);
    printf("%s = %s\n", tr.getTag()[0].c_str(), tr.getAnnotation()[0].c_str());
    printf("%s = %s\n", tr.getTag()[1].c_str(), tr.getAnnotation()[1].c_str());
  }

  {
    const char path[] = "/net/fantasia/home/hmkang/bin/annovar/humandb/hg19_ljb_all.txt.gz";
    TabixReader tr(path, 1, 2, 4, 5); // dbNSFP chromosome
    
    tr.addTag("SIFT", 6);
    tr.addTag("GERP", 16);
    tr.addAnnotation("1", 69091, "A", "C");
    assert(tr.getAnnotation().size() == 2);
    printf("%s = %s\n", tr.getTag()[0].c_str(), tr.getAnnotation()[0].c_str());
    printf("%s = %s\n", tr.getTag()[1].c_str(), tr.getAnnotation()[1].c_str());
  }
  
  return 0;
}
