all: Main
Main: Main.cpp 
	g++ -g -c Main.cpp -I/home/zhanxw/statgen/lib/general
	g++ -g -o Main Main.o ../statgen/lib/libStatGen.a -lz -lbz2 -lssl
clean:
	rm -f *.o Main
