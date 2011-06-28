EXEC = Main
DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS

.PHONY: release debug
all: debug
release: CXXFLAGS = -O2 $(DEFAULT_CXXFLAGS)
release: $(EXEC)
debug: CXXFLAGS = -g $(DEFAULT_CXXFLAGS)
debug: $(EXEC)

Main: Main.cpp Gene.h Range.h IO.h Argument.h FreqTable.h
	g++ $(CXXFLAGS) -c Main.cpp
	g++ $(CXXFLAGS) -o Main Main.o -lz -lbz2 -lssl -lcrypto
clean:
	rm -f *.o Main
test1: debug
	./Main -i test.vcf -r test.fa -g test.gene.txt -o test.out.vcf
test2: debug
	./Main -i 100.vcf.gz -r test.fa -g refFlat_hg19.txt.gz -o 100.out.vcf
Log: LogFile.cpp LogFile.h
	g++ -g -o $@ $<
