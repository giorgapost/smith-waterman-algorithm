@echo off

echo Building object files...
IF not exist obj (md obj)
g++ -Wall -fexceptions -g -Isrc\include -c src\SmithWatermanExecutor.cpp -o obj\SmithWatermanExecutor.o
g++ -Wall -fexceptions -g -Isrc\include -c src\classes\Framework.cpp -o obj\Framework.o
g++ -Wall -fexceptions -g -Isrc\include -c src\classes\ParallelCoarseOMPImplementation.cpp -o obj\ParallelCoarseOMPImplementation.o -fopenmp
g++ -Wall -fexceptions -g -Isrc\include -c src\classes\ParallelFineOMPImplementation.cpp -o obj\ParallelFineOMPImplementation.o -fopenmp
g++ -Wall -fexceptions -g -Isrc\include -c src\classes\SequentialImplementation.cpp -o obj\SequentialImplementation.o

echo Linking files and constructing an executable...
g++  -o smith_waterman.exe obj\SmithWatermanExecutor.o obj\Framework.o obj\ParallelCoarseOMPImplementation.o obj\ParallelFineOMPImplementation.o obj\SequentialImplementation.o -fopenmp
IF not exist reports (md reports)

echo Cleaning up...
rmdir /s /Q obj
