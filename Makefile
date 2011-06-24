EXEC = Main
DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS

.PHONY: release debug
all: debug
release: CXXFLAGS = -O2 $(DEFAULT_CXXFLAGS)
release: $(EXEC)
debug: CXXFLAGS = -g $(DEFAULT_CXXFLAGS)
debug: $(EXEC)

Main: Main.cpp Gene.h Range.h
	g++ $(CXXFLAGS) -c Main.cpp -I../statgen/lib/general
	g++ $(CXXFLAGS) -o Main Main.o ../statgen/lib/libStatGen.a -lz -lbz2 -lssl -lcrypto
clean:
	rm -f *.o Main
test1: Main
	./Main -i test.vcf -r test.fa -g test.gene.txt
test2: Main
	./Main -i 100.vcf -r test.fa -g refFlat_hg19.txt.gz
