/**
 * @file Framework.h
 */
#pragma once

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <sys/time.h>

using namespace std;

/**
* A pair of 2 sequences named Q and D, which have to be aligned
* (by the <a href="https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm">Smith-Waterman algorithm</a>).
* @author Georgios Apostolakis
*/
struct Pair {
    /** A string with the Q sequence. */
    string q;

    /** A string with the D sequence. */
    string d;
};

/**
* The coordinates of a 2D matrix cell.
* @author Georgios Apostolakis
*/
struct Position {
    /** The row of the cell. */
    unsigned int row;

    /** The column of the cell. */
    unsigned int column;
};

/**
* A possible alignment for an input pair of sequences.
* @author Georgios Apostolakis
*/
struct Result {
    /** An integer ID that indicates the pair of sequences to which this object corresponds. */
    unsigned int ref_id;

    /** The score of the current result. */
    int score;

    /** An integer with the starting index. */
    int start;

    /** An integer with the stopping index. */
    int stop;

    /** A {@link #Pair} object with the aligned sequences. */
    Pair result_pair;
};

/**
* The parameters that define the behaviour and output of the algorithm.
* @author Georgios Apostolakis
*/
struct Scores {
    /** The score of a match between the sequences. */
    int matchScore;

    /** The score of a mismatch between the sequences. */
    int mismatchScore;

    /** The score of a gap between the sequences. */
    int gapScore;
};

/**
* Some statistical data about the algorithm's execution.
* @author Georgios Apostolakis
*/
struct Statistics {
    /** The cells at the scoring matrix with a value greater than zero. */
    long long int cellsGreaterThanZero;

    /** The total traceback steps executed by the algorithm. */
    long long int totalTracebackSteps;

    /** The total time that the execution of the algorithm lasted (in seconds). */
    long double totalTime;

    /** The total time that the computation of the scoring matrix cells lasted (in seconds). */
    long double *calcCellsTime;

    /** The total time that the traceback process lasted (in seconds). */
    long double *totalTracebackTime;
};

/**
 * This abstract class can be used as a framework which provides functions that facilitate
 * various implementations of the
 * <a href="https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm">Smith-Waterman algorithm</a>.
 * @author Georgios Apostolakis
 */
class Framework {
    public:
        /**
         * Constructs a new instance of this class. If no arguments are provided by the user,
         * it asks for them through console questions. Then, it reads the sequences to be aligned
         * from an input file.<br>
         * Valid arguments:
         * - {@code -id <string with the id of the report>}
         * - {@code -path <string with the input file's path>}
         * - {@code -match <integer with the match score>}
         * - {@code -mismatch <integer with the mismatch score>}
         * - {@code -gap <integer with the gap score>}<br>
         * Notice that any extra (and possibly invalid) arguments are ignored without throwing any exception.
         * @param argc An integer with the size of the {@code argv} argument.
         * @param argv An array with the arguments provided by the user.
         * @throws std::ios_base::failure Thrown if the input file cannot be opened.
         * @throws std::invalid_argument Thrown if some arguments from the listed above are missing.
         * @throws std::runtime_error Thrown if the contents of the input file are invalid.
         */
        Framework(int argc, char* argv[]);

        /**
         * Destroys an instance of the current class.
         */
        virtual ~Framework(void);

        /**
         * Saves the aligned sequences into file '/reports/Report_ID.txt',
         * where ID was given by the user (either as an argument, or through console).
         * The data to save is retrieved from the {@link Framework#results results} member-variable of this class.
         * @throws std::ios_base::failure Thrown if the output file cannot be opened.
         */
        void printResultsToFile(void);

        /**
         * An abstract method that prints some statistics into console about the execution details of the algorithm.
         */
        virtual void printStatistics(void) = 0;

        /**
         * This abstract method executes the
         * <a href="https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm">Smith-Waterman algorithm</a>
         * and has to be implemented in all descendant classes.
         */
        virtual void runAlgorithm(void) = 0;

    protected:
        /** Contains the match, mismatch, gap scores on which depends the output of the algorithm. */
        Scores algoScores;

        /** A vector of {@link Pair} objects, with pairs of sequences to be aligned. */
        vector<Pair> data;

        /** The path of the file that contains the input data, i.e. the sequences that need alignment. */
        string path;

        /**
         * A string with the ID of the report produced after the execution of the
         * algorithm. The ID is part of the report's filename.
         */
        string reportId;

        /** An "array" of {@link Result} objects, with all the optimal alignments for every input pair of sequences. */
        vector<vector<Result>> results;

        /** Contains statistical data about the execution of the algorithm. */
        Statistics statisticData;

        /**
         * Provides the current time (since the
         * <a href="https://en.wikipedia.org/wiki/Epoch_(computing)">Epoch</a>) in seconds.
         * @return A long double variable with the number of seconds elapsed since the Epoch.
         */
        long double getTime(void);

    private:
        /**
         * Reads the required arguments from console, by making appropriate questions to the user.
         */
        void readArgsFromConsole(void);

        /**
         * Reads the input sequences (which need alignment) from a file whose path is given by
         * {@link #path} member-variable.
         * It stores them into the {@link #data} member-variable.
         */
        void readInputFile(void);
};
