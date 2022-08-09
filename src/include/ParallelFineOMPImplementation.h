/**
 * @file ParallelFineOMPImplementation.h
 */
#pragma once

#include "Framework.h"
#include <omp.h>

using namespace std;

/**
 * This class extends the {@link Framework} class and implements
 * the <a href="https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm">Smith-Waterman algorithm</a>
 * in a parallel multi-threaded way. Parallelization takes place in a fine level, where the cells of
 * the scoring matrix for every pair of sequences are computed in parallel.
 * @author Georgios Apostolakis
 */
class ParallelFineOMPImplementation : public Framework{
    public:
        /**
         * Calls the Framework() constructor of the parent class and also determines the number of threads that will be used.
         * If no arguments are provided by the user it asks for them through console questions.
         * The valid arguments are the same with the arguments listed in the documentation of the Framework() constructor, plus:
         * - {@code -threads <integer with the threads for parallelization>}<br>
         * Any extra (and possibly invalid) arguments are ignored without throwing any exception.
         * @param argc An integer with the size of the {@code argv} argument.
         * @param argv An array with the arguments provided by the user.
         * @throws std::ios_base::failure Thrown if the input or output file cannot be opened.
         * @throws std::invalid_argument Thrown if some arguments from the listed above are missing.
         * @throws std::runtime_error Thrown if the contents of the input file are invalid.
         */
        ParallelFineOMPImplementation(int argc, char* argv[]);

        /**
         * Destroys an instance of the current class by calling the ~Framework() destructor
         * of the parent class.
         */
        virtual ~ParallelFineOMPImplementation();

        /**
         * Executes the <a href="https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm">Smith-Waterman algorithm</a>
         * in a parallel way (multiple threads). Parallelization takes place through the computation of different cells
         * of the scoring matrix. However, all the alignments for a given pair have to be computed before the implementation
         * moves on to the next one (on the contrary to the coarse level parallelization).
         */
        void runAlgorithm(void) final;

        /**
         * Prints some statistics into console about the execution details of the algorithm.
         */
        void printStatistics(void) final;

    private:
        /** The maximum number of threads to use for parallelization **/
        unsigned int threads;

        /**
         * Fills the scoring matrix of the
         * <a href="https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm">Smith-Waterman algorithm</a>
         * in a parallel multi-threaded way (for the {@link Pair} of sequences provided as argument). Different
         * cells may be computed by different threads.
         * @param sequences A {@link Pair} object with the sequences that will be aligned.
         * @param scoring_matrix A 2D matrix constructed of vectors, whose cells will be filled according to the process
         * defined by the algorithm. Every cell will contain an integer value.
         * @param max_pos_vec A vector whose entries are of type {@link Position} and hold the coordinates of the cells with maximum value.
         * @return A long long integer with the number of cells whose entries are greater than zero.
         */
        long long int fill_scoring_matrix(Pair sequences, vector< vector<int> >& scoring_matrix, vector<Position>& max_pos_vec);

        /**
         * Performs a recursive process in order to extract the aligned sequences from the scoring matrix.
         * @param steps The number of traceback steps performed until the beginning of the recursive method, increased by 1.
         * @param start_row The row of the cell where the traceback will begin.
         * @param start_col The column of the cell where the traceback will begin.
         * @param sequences A {@link Pair} object with the sequences that will be aligned.
         * @param scoring_matrix The scoring matrix constructed by the #fill_scoring_matrix() method.
         * @param results_vec A vector whose entries are of type {@link Result} and holds all the optimal alignments for the given sequences.
         * @return A long long integer with the number of total traceback steps that were required (i.e. the depth of the recursion).
         */
        long long int traceback(int steps, int start_row, int start_col, Pair sequences, vector< vector<int> >& scoring_matrix, vector<Result>& results_vec);
};
