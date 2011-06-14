all: Main
Main: Main.cpp
	g++ -o $@ $<
