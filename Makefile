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

$(EXEC): Main.cpp Gene.h Range.h IO.h Argument.h FreqTable.h GenomeSequence.h
	g++ $(CXXFLAGS) -c Main.cpp
	g++ $(CXXFLAGS) -o $@ Main.o -lz -lbz2 -lssl -lcrypto
clean:
	rm -f *.o $(EXEC)
# basic test
test1: debug
	(cd example; ../$(EXEC) -i test.vcf -r test.fa -g test.gene.txt -c ../codon.txt -o test.out.vcf)
test2: debug
	(cd example; ../$(EXEC) -i 100.vcf.gz -r test.fa -g ../resources/refFlat_hg19.txt.gz -c ../codon.txt -o noref.100.out.vcf)
test3: debug
	(cd example; ../$(EXEC) -i 100.vcf.gz -r ../resources/human.g1k.v37.fa -g ../resources/refFlat_hg19.txt.gz -c ../codon.txt -o 100.out.vcf)

# test plain text file annotation
test4: debug
	(cd example; ../$(EXEC) -i test.plain.txt -r test.fa -g test.gene.txt -c ../codon.txt -o test.plain.anno.txt --inputFormat plain)

test: test1 test2 test3 test4
check:
	diff example/test.out.vcf example/correct.test.out.vcf && \
	diff example/noref.100.out.vcf example/correct.noref.100.out.vcf && \
	diff example/100.out.vcf example/correct.100.out.vcf && \
	diff example/test.plain.anno.txt example/correct.test.plain.anno.txt 

correct:
	@echo "regenerate correct testing files..."
	sleep 10
	(cd example; ../executable/$(EXEC) -i test.vcf -r test.fa -g test.gene.txt -c ../codon.txt -o correct.test.out.vcf)
	(cd example; ../executable/$(EXEC) -i 100.vcf.gz -r test.fa -g ../resources/refFlat_hg19.txt.gz -c ../codon.txt -o correct.noref.100.out.vcf)
	(cd example; ../executable/$(EXEC) -i 100.vcf.gz -r ../resources/human.g1k.v37.fa -g ../resources/refFlat_hg19.txt.gz -c ../codon.txt -o correct.100.out.vcf)
	(cd example; ../executable/$(EXEC) -i test.plain.txt -r test.fa -g test.gene.txt -c ../codon.txt -o correct.test.plain.anno.txt)

# auxillary tools
Log: LogFile.cpp LogFile.h
	g++ -g -o $@ $<
