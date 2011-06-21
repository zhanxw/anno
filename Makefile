all: Main
Main: Main.cpp 
	g++ -c Main.cpp -I/home/zhanxw/statgen/lib/general
	g++ -o Main Main.o ../statgen/lib/libStatGen.a -lz -lbz2 -lssl
