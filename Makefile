CXXFLAGS = -g
all: Main
Main: Main.cpp Gene.h Range.h
	g++ $(CXXFLAGS) -c Main.cpp -I../statgen/lib/general
	g++ $(CXXFLAGS) -o Main Main.o ../statgen/lib/libStatGen.a -lz -lbz2 -lssl -lcrypto
clean:
	rm -f *.o Main
