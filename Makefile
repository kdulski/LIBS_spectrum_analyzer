CC = g++

CFLAGS = -Wall `root-config --cflags`

CLIBS = `root-config --glibs` `root-config --cflags`

all: LIBSanalyzer

LIBSanalyzer: spectrum_analysis.cpp
	$(CC) spectrum_analysis.cpp -o LIBSAna.exe $(CLIBS) -std=c++17

clean: 
	rm -rf *o

