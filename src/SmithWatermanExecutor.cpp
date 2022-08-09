/**
 * @file SmithWatermanExecutor.cpp
 */
#include "SequentialImplementation.h"
#include "ParallelCoarseOMPImplementation.h"
#include "ParallelFineOMPImplementation.h"

using namespace std;

/** The integer value of this constant is the value of the input argument '-parallel' that determines a sequential execution of the algorithm.*/
const int SEQUENTIAL_IMPL = 1;

/** The integer value of this constant is the value of the input argument '-parallel' that determines a coarse-level parallel execution of the algorithm.*/
const int PARALLEL_COARSE_IMPL = 2;

/** The integer value of this constant is the value of the input argument '-parallel' that determines a fine-level parallel execution of the algorithm.*/
const int PARALLEL_FINE_IMPL = 3;

/**
 * Determines which version of the <a href="https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm">Smith-Waterman algorithm</a> to execute, based
 * on the input arguments. If no arguments were provided, it lets the user to
 * choose an option through console.
 * @param argc An integer with the size of the {@code argv} argument.
 * @param argv An array with the arguments provided by the user.
 * @return An integer with the numerical value of the respective constant
 * ({@link #SEQUENTIAL_IMPL}, {@link #PARALLEL_COARSE_IMPL}, {@link #PARALLEL_FINE_IMPL})
 * which indicates the version of the algorithm to execute.
 */
int selectAlgorithm(int argc, char* argv[]){
    int algo = -1;
    if(argc<=1){
        cin.clear();
        cout << "Please enter the number of the implementation to execute:" << endl;
        cout << "  " + to_string(SEQUENTIAL_IMPL) + ". Sequential implementation." << endl;
        cout << "  " + to_string(PARALLEL_COARSE_IMPL) + ". Parallel coarse-grained implementation." << endl;
        cout << "  " + to_string(PARALLEL_FINE_IMPL) + ". Parallel fine-grained implementation." << endl;
        cin >> algo;
    }
    else {
        for(int i=0;i<argc-1;i++){
            if(!string(argv[i]).compare("-parallel")){
                algo = atoi(argv[i+1]);
                break;
            }
        }
    }
    return algo;
}

/**
 * The main function, which initializes the execution of the appropriate version of the <a href="https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm">Smith-Waterman algorithm</a>.
 * @param argc An integer with the size of the {@code argv} argument.
 * @param argv An array with the arguments provided by the user.
 * @return An integer with the value 0.
 */
int main(int argc, char* argv[]){
    int algo = selectAlgorithm(argc, argv);
    if(algo==SEQUENTIAL_IMPL){
        try{
            SequentialImplementation ser(argc, argv);
            ser.runAlgorithm();
            ser.printResultsToFile();
            ser.printStatistics();
        }catch(const std::exception& e) {
            cerr << e.what() << " Program will be terminated." << endl;
        }
    }
    else if(algo==PARALLEL_COARSE_IMPL){
        try{
            ParallelCoarseOMPImplementation par(argc, argv);
            par.runAlgorithm();
            par.printResultsToFile();
            par.printStatistics();
        }catch(const std::exception& e) {
            cerr << e.what() << " Program will be terminated." << endl;
        }
    }
    else if(algo==PARALLEL_FINE_IMPL){
        try{
            ParallelFineOMPImplementation par(argc, argv);
            par.runAlgorithm();
            par.printResultsToFile();
            par.printStatistics();
        }catch(const std::exception& e) {
            cerr << e.what() << " Program will be terminated." << endl;
        }
    }
    else
		cerr << "Error. Invalid arguments. Program will be terminated." << endl;
	
	return 0;
}
