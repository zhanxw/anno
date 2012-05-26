EXEC = anno
DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS

.PHONY: release debug
all: debug
release: CXXFLAGS = -O2 $(DEFAULT_CXXFLAGS)
release: $(EXEC)
	-mkdir -p executable
	cp -f $(EXEC) executable
debug: CXXFLAGS = -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(EXEC)

$(EXEC): Main.cpp Gene.h Range.h IO.h Argument.h FreqTable.h GenomeSequence.h LogFile.h StringTemplate.h
	g++ $(CXXFLAGS) -c Main.cpp
	g++ $(CXXFLAGS) -o $@ Main.o -lz -lbz2 -lssl -lcrypto
clean:
	rm -f *.o $(EXEC) input.*
# basic test
test1: debug
	(cd example; ../$(EXEC) -i input.test.vcf -r test.fa -g test.gene.txt -c ../codon.txt -o output.test.vcf)
	(cd example; diff {correct,output}.test.vcf)

test2: debug
	(cd example; ../$(EXEC) -i input.100.vcf.gz -r test.fa -g ../resources/refFlat_hg19.txt.gz -c ../codon.txt -o output.noref.100.vcf)
	(cd example; diff {correct,output}.noref.100.vcf)

test3: debug
	(cd example; ../$(EXEC) -i input.100.vcf.gz -r ../resources/human.g1k.v37.fa -g ../resources/refFlat_hg19.txt.gz -c ../codon.txt -o output.100.vcf)
	(cd example; diff {correct,output}.100.vcf)

# test plain text file annotation
test4: debug
	(cd example; ../$(EXEC) --inputFormat plain -i input.test.plain.txt -r test.fa -g test.gene.txt -c ../codon.txt -o output.test.plain.txt --inputFormat plain)
	(cd example; diff {correct,output}.test.plain.txt)

# test refGene format
test5: debug
	(cd example; ../$(EXEC) -i input.100.vcf.gz -r ../resources/human.g1k.v37.fa -f refGene -g ../resources/refGene.txt.gz -c ../codon.txt -o output.100.refGene.vcf)
	(cd example; diff {correct,output}.100.refGene.vcf)

# test epacts annotation
test6: debug
	(cd example; ../$(EXEC) -i input.epacts --inputFormat epacts -r ../resources/human.g1k.v37.fa -f refGene -g ../resources/refGene.txt.gz -c ../codon.txt -o output.epacts)
	(cd example; diff {correct,output}.epacts)

test: test1 test2 test3 test4 test5 test6

correct:
	@echo "WARNING: Regenerate correct testing files. You have 10 seconds to stop..."
	sleep 10
	(cd example; ../$(EXEC) -i input.test.vcf -r test.fa -g test.gene.txt -c ../codon.txt -o correct.test.vcf)
	(cd example; ../$(EXEC) -i input.100.vcf.gz -r test.fa -g ../resources/refFlat_hg19.txt.gz -c ../codon.txt -o correct.noref.100.vcf)
	(cd example; ../$(EXEC) -i input.100.vcf.gz -r ../resources/human.g1k.v37.fa -g ../resources/refFlat_hg19.txt.gz -c ../codon.txt -o correct.100.vcf)
	(cd example; ../$(EXEC) --inputFormat plain -i input.test.plain.txt -r test.fa -g test.gene.txt -c ../codon.txt -o correct.test.plain.txt --inputFormat plain)
	(cd example; ../$(EXEC) -i input.100.vcf.gz -r ../resources/human.g1k.v37.fa -f refGene -g ../resources/refGene.txt.gz -c ../codon.txt -o correct.100.refGene.vcf)
	(cd example; ../$(EXEC) -i input.epacts --inputFormat epacts -r ../resources/human.g1k.v37.fa -f refGene -g ../resources/refGene.txt.gz -c ../codon.txt -o correct.epacts)
	# (cd example; ../executable/$(EXEC) -i input.test.vcf -r test.fa -g test.gene.txt -c ../codon.txt -o correct.test.vcf)
	# (cd example; ../executable/$(EXEC) -i input.100.vcf.gz -r test.fa -g ../resources/refFlat_hg19.txt.gz -c ../codon.txt -o correct.noref.100.vcf)
	# (cd example; ../executable/$(EXEC) -i input.100.vcf.gz -r ../resources/human.g1k.v37.fa -g ../resources/refFlat_hg19.txt.gz -c ../codon.txt -o correct.100.vcf)
	# (cd example; ../executable/$(EXEC) --inputFormat plain -i input.test.plain.txt -r test.fa -g test.gene.txt -c ../codon.txt -o correct.test.plain.txt)


# auxillary tools
Log: LogFile.cpp LogFile.h
	g++ -g -o $@ $<


testStringTemplate: testStringTemplate.cpp StringTemplate.h
	g++ -g -o $@ $<

doc:
	java -jar ext/wiki2html.jar README.wiki > README.html 
	pandoc -f html -t markdown README.html > README.md 
