#GWAMA program

VERSION = 2.7

CC = g++

DEBUGFLAGS = -Wno-deprecated -O3 -lz

GWAMA:	main.cpp

	g++ main.cpp marker.cpp statistics.cpp study.cpp chisquaredistr.cpp normaldistr.cpp gammaf.cpp igammaf.cpp ap.cpp global.cpp problem.cpp tools.cpp cohort.cpp commandLine.cpp readFile.cpp $(DEBUGFLAGS) -o GWAMA

