EXEC = anno
DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS

.PHONY: release debug
all: debug
release: CXXFLAGS = -O2 $(DEFAULT_CXXFLAGS)
release: $(EXEC)
	-mkdir executable
	cp -f $(EXEC) executable
debug: CXXFLAGS = -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(EXEC)

$(EXEC): Main.cpp Gene.h Range.h IO.h Argument.h FreqTable.h
	g++ $(CXXFLAGS) -c Main.cpp
	g++ $(CXXFLAGS) -o $@ Main.o -lz -lbz2 -lssl -lcrypto
clean:
	rm -f *.o $(EXEC)
test1: debug
	(cd example; ../$(EXEC) -i test.vcf -r test.fa -g test.gene.txt -c ../codon.txt -o test.out.vcf)
test2: debug
	(cd example; ../$(EXEC) -i 100.vcf.gz -r test.fa -g refFlat_hg19.txt.gz -c ../codon.txt -o 100.out.vcf)
Log: LogFile.cpp LogFile.h
	g++ -g -o $@ $<
