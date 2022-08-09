# the compiler: define as g++ for C++
CXX = g++

# Definition of paths
OBJECT_DIR = obj/
REPORT_DIR = reports/
INCLUDE_DIR = src/include
SOURCE_DIR = src/
CLASSES_DIR = src/classes/

# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS  = -g -Wall -fexceptions -I$(INCLUDE_DIR)

# Necessary libraries
LIBS = -fopenmp

# The build target 
TARGET = smith_waterman

.PHONY: all clean $(TARGET) obj_files
 
all: $(TARGET)
 
$(TARGET): obj_files
	@mkdir -p $(OBJECT_DIR)
	$(CXX) -o $(TARGET) $(OBJECT_DIR)SmithWatermanExecutor.o $(OBJECT_DIR)Framework.o $(OBJECT_DIR)ParallelCoarseOMPImplementation.o $(OBJECT_DIR)ParallelFineOMPImplementation.o $(OBJECT_DIR)SequentialImplementation.o $(LIBS)
	@mkdir -p $(REPORT_DIR)
		
obj_files:
	@mkdir -p $(OBJECT_DIR)
	$(CXX) $(CFLAGS) -c $(SOURCE_DIR)SmithWatermanExecutor.cpp -o $(OBJECT_DIR)SmithWatermanExecutor.o
	$(CXX) $(CFLAGS) -c $(CLASSES_DIR)Framework.cpp -o $(OBJECT_DIR)Framework.o
	$(CXX) $(CFLAGS) -c $(CLASSES_DIR)ParallelCoarseOMPImplementation.cpp -o $(OBJECT_DIR)ParallelCoarseOMPImplementation.o $(LIBS)
	$(CXX) $(CFLAGS) -c $(CLASSES_DIR)ParallelFineOMPImplementation.cpp -o $(OBJECT_DIR)ParallelFineOMPImplementation.o $(LIBS)
	$(CXX) $(CFLAGS) -c $(CLASSES_DIR)SequentialImplementation.cpp -o $(OBJECT_DIR)SequentialImplementation.o $(LIBS)

clean:
	rm -rf $(OBJECT_DIR)
	rm -f $(TARGET)