EXEC = anno
DEFAULT_CXXFLAGS = -D__STDC_LIMIT_MACROS


.PHONY: release debug
all: debug

LIB = third/tabix/libtabix.a
third/tabix/libtabix.a:
	(cd third; make tabix)

release: CXXFLAGS = -O2 -DNDEBUG $(DEFAULT_CXXFLAGS)
release: $(EXEC)
	-mkdir -p executable
	cp -f $(EXEC) executable
debug: CXXFLAGS = -Wall -ggdb -O0 $(DEFAULT_CXXFLAGS)
debug: $(EXEC)
profile: CXXFLAGS = -ggdb -pg -O0 $(DEFAULT_CXXFLAGS)
profile: $(EXEC)

GitVersion.h: .git/HEAD .git/index
	echo "const char *gitVersion = \"$(shell git rev-parse HEAD)\";" > $@

TabixReader.h: $(LIB)

-include IO.d
IO.o: IO.h IO.cpp
	g++ -MMD $(CXXFLAGS) -c IO.cpp
-include Main.d
Main.o: Main.cpp $(LIB)
	g++ -MMD $(CXXFLAGS) -c Main.cpp

$(EXEC): Main.o IO.o
	g++ $(CXXFLAGS) -o $@ Main.o IO.o $(LIB) -lz -lbz2 -lssl -lcrypto
clean:
	rm -f *.o *.d $(EXEC) input.*

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
	(cd example; ../$(EXEC) -i input.epacts --inputFormat epacts -r ../resources/human.g1k.v37.fa -f refGene -g ../resources/refGene.txt.gz -c ../codon.txt -o output.epacts --outputFormat epact -p ../priority.epact.txt)
	(cd example; diff {correct,output}.epacts)

test7: debug
	(cd example; ../$(EXEC) -i input.test.vcf -r test.fa -g test.gene.txt -c ../codon.txt -o output.test.region.vcf --bed REGION1=region1.bed)
#	(cd example; diff {correct,output}.test.vcf)

test8: debug
	(cd example; ../$(EXEC) -i input.test.vcf -r test.fa -g test.gene.txt -c ../codon.txt -o output.test.region.vcf --bed REGION1=region1.bed,REGION2=region2.bed)

testBedReader: testBedReader.cpp BedReader.h
	g++ $(CXXFLAGS) -c testBedReader.cpp
	g++ $(CXXFLAGS) -o $@ testBedReader.o -lz -lbz2 -lssl -lcrypto

test9: testBedReader
	(cd example; ../testBedReader > output.testBedReader)
	(cd example; diff  {correct,output}.testBedReader)

testGenomeScore: testGenomeScore.cpp GenomeScore.h
	g++ $(CXXFLAGS) -c testGenomeScore.cpp -O2
	g++ $(CXXFLAGS) -o $@ testGenomeScore.o -lz -lbz2 -lssl -lcrypto -O2

test10: testGenomeScore
	(cd example; ../testGenomeScore > output.testGenomeScore)
	(cd example; diff {correct,output}.testGenomeScore)

testTabixReader: testTabixReader.cpp TabixReader.h
	g++ $(CXXFLAGS) -c testTabixReader.cpp -O0 -ggdb
	g++ $(CXXFLAGS) -o $@ testTabixReader.o $(LIB) -lz -lbz2 -lssl -lcrypto

test11: testTabixReader 
	(cd example; ../testTabixReader > output.testTabixReader)
	(cd example; diff {correct,output}.testTabixReader)

# test gz output and index
test12: debug
	(cd example; ../$(EXEC) -i input.100.vcf.gz -r ../resources/human.g1k.v37.fa -f refGene -g ../resources/refGene.txt.gz -c ../codon.txt -o output.100.refGene.vcf.gz --indexOutput)
	(cd example; diff {correct,output}.100.refGene.vcf.gz; diff {correct,output}.100.refGene.vcf.gz.tbi)

# test gz output for plain text file annotation
test13: debug
	(cd example; ../$(EXEC) --inputFormat plain -i input.test2.plain.txt -r test.fa -g test.gene.txt -c ../codon.txt -o output.test2.plain.txt.gz --indexOutput)
	(cd example; diff {correct,output}.test2.plain.txt.gz; diff {correct,output}.test2.plain.txt.gz.tbi)

# test gz output and index for plain text file annotation
test14: debug
	(cd example; ../$(EXEC) --inputFormat plain -i input.test3.plain.txt -r test.fa -g test.gene.txt -c ../codon.txt -o output.test3.plain.txt.gz --indexOutput)
	(cd example; zdiff {correct,output}.test3.plain.txt.gz; diff {correct,output}.test3.plain.txt.gz.tbi)

test15: testLineBreaker

test: test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test11 test12 test13 test14 test15


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
testLineBreaker: testLineBreaker.cpp LineBreaker.h
	g++ -g -o $@ $<
doc:
	java -jar ext/wiki2html.jar README.wiki > README.html 
	pandoc -f html -t markdown README.html > README.md 
