/**
 * @file ParallelFineOMPImplementation.cpp
 */
#include "ParallelFineOMPImplementation.h"

ParallelFineOMPImplementation::ParallelFineOMPImplementation(int argc, char* argv[]): Framework(argc, argv) {
    int i;
    if(argc<=1){
        cin.clear();
        cout << "Please enter the number of threads for parallel execution:";
        cin >> threads;
    }
    else{
         for(i=0;i<argc-1;i+=1){
            if(!string(argv[i]).compare("-threads")){
                threads = atoi(argv[i+1]);
                break;
            }
         }
         if(i==argc-1)
            throw invalid_argument("Error. Missing arguments.");
    }

    statisticData.calcCellsTime = new long double[1];
    statisticData.totalTracebackTime = new long double[1];
    statisticData.calcCellsTime[0] = 0;
    statisticData.totalTracebackTime[0] = 0;
}

ParallelFineOMPImplementation::~ParallelFineOMPImplementation(){
    delete statisticData.calcCellsTime;
    delete statisticData.totalTracebackTime;
}

long long int ParallelFineOMPImplementation::fill_scoring_matrix(Pair sequences, vector< vector<int> >& scoring_matrix, vector<Position>& max_pos_vec){
	vector<Position> max_vecs[threads];
	int max_vals[threads], cell_vals[threads];;

    #pragma omp parallel shared(sequences, scoring_matrix, statisticData, max_vecs, max_vals, cell_vals) num_threads(threads)
    {
        unsigned int row = 0, col = 0;
        int size0, i, cellCounter=0, n1, n2, n3;
        int maxNum=0;
        vector<Position> max_pos;

        unsigned int num_rows = sequences.q.size()+1;  //number of rows of the matrix
        unsigned int num_cols = sequences.d.size()+1;  //number of columns of the matrix

        while(row<num_rows && col<num_cols){  //max value of row will be num_rows-1, but max value of col will be num_cols
            //calculate in parallel all cells for the minor diagon of the matrix which starts at (row, 0) or (max_row, col)

            if(row<num_rows-1)
                size0 = row+1;
            else if(num_cols-col<=num_rows)
                size0 = num_cols-col;
            else
                size0 = min(num_rows, num_cols);

            #pragma omp for
            for (i=0; i<size0; i++) {

                //calculate scoring_matrix[row-i][col-i]
                if((row-i==0) || (col+i==0))
                    scoring_matrix[row-i][col+i] = 0;
                else {
                    if(sequences.q[row-i-1]==sequences.d[col+i-1])
                        n1 = scoring_matrix[row-i-1][col+i-1] + algoScores.matchScore;
                    else
                        n1 = scoring_matrix[row-i-1][col+i-1] + algoScores.mismatchScore;

                    n2 = scoring_matrix[row-i-1][col+i] + algoScores.gapScore; //up
                    n3 = scoring_matrix[row-i][col+i-1] + algoScores.gapScore; //left

                    scoring_matrix[row-i][col+i] = max(max(0, n1), max(n2, n3));
                }

                Position tmp_pos;
                if(scoring_matrix[row-i][col+i]>maxNum){
                    max_pos.clear();
                    tmp_pos.row=row-i;
                    tmp_pos.column=col+i;
                    max_pos.push_back(tmp_pos);
                    maxNum = scoring_matrix[row-i][col+i];
                }
                else if(scoring_matrix[row-i][col+i]==maxNum){
                    tmp_pos.row=row-i;
                    tmp_pos.column=col+i;
                    max_pos.push_back(tmp_pos);
                }

                if(scoring_matrix[row-i][col+i]>0)
                    cellCounter++;
            }

            if(row<num_rows-1)
                row++;
            else
                col++;
        }

        max_vecs[omp_get_thread_num()] = max_pos;
        max_vals[omp_get_thread_num()] = maxNum;
        cell_vals[omp_get_thread_num()] = cellCounter;
    }
	
    int max_v=max_vals[0];
	for(unsigned int i=0;i<threads;i++){
		if(max_vals[i]>max_v){
			max_v = max_vals[i];
			max_pos_vec = max_vecs[i];
		}
		else if(max_vals[i]==max_v)
			for(unsigned int j=0;j<max_vecs[i].size();j++)
				max_pos_vec.push_back(max_vecs[i][j]);
	}

	for(unsigned int i=0;i<threads;i++)
		statisticData.cellsGreaterThanZero += cell_vals[i];

	return statisticData.cellsGreaterThanZero;
}

void ParallelFineOMPImplementation::printStatistics(void){
    cout << "A) Total pairs of sequences Q-D: " << data.size() << endl;
    cout << "B) Total cells with value: " << statisticData.cellsGreaterThanZero << endl;
    cout << "C) Total traceback steps: " << statisticData.totalTracebackSteps << endl;
    cout << "D) Total time of program execution: " << statisticData.totalTime << " seconds"<< endl;
    cout << "E) Total time of calculating cells: " << *statisticData.calcCellsTime << " seconds" << endl;
    cout << "F) Total traceback time: " << *statisticData.totalTracebackTime << " seconds" << endl;
    cout << "G) Cell Updates Per Second (CUPS) based on total execution time: " << (double)statisticData.cellsGreaterThanZero/statisticData.totalTime << endl;
    cout << "H) Cell Updates Per Second (CUPS) based on time of calculating cells: " << (double)statisticData.cellsGreaterThanZero/(*statisticData.calcCellsTime) << endl;
}

void ParallelFineOMPImplementation::runAlgorithm(void){
    long double time_t0 = getTime();

    for(size_t i=0;i<data.size();i++){
        vector< vector<int> > scoring_matrix(data[i].q.size()+1, vector<int>(data[i].d.size()+1));
        vector<Position> max_pos_vec;

		long double time_f1 = getTime();  //Filling the scoring matrix
		statisticData.cellsGreaterThanZero += fill_scoring_matrix(data[i], scoring_matrix, max_pos_vec);
        *statisticData.calcCellsTime += getTime() - time_f1;


        for(size_t j=0;j<max_pos_vec.size();j++){
            Result res;
            res.ref_id = i;
			res.score = scoring_matrix[max_pos_vec[j].row][max_pos_vec[j].column];
			res.stop = max_pos_vec[j].column-1;
			results[i].push_back(res);

            long double time_tr1 = getTime();
			statisticData.totalTracebackSteps += traceback(1, max_pos_vec[j].row, max_pos_vec[j].column, data[i], scoring_matrix, results[i]);
			*statisticData.totalTracebackTime += getTime() - time_tr1;
		}
	}
	statisticData.totalTime += getTime() - time_t0;
}
long long int ParallelFineOMPImplementation::traceback(int steps, int start_row, int start_col, Pair sequences, vector< vector<int> >& scoring_matrix, vector<Result>& results_vec){
	unsigned int index = results_vec.size()-1;
	int n1, n2, n3;

	//Compute values of left, up, and diagonally left cell of cell [start_row,start_col]
	if(sequences.q[start_row-1]==sequences.d[start_col-1])			 //diagonally up left
		n1 = scoring_matrix[start_row-1][start_col-1] + algoScores.matchScore;
	else
		n1 = scoring_matrix[start_row-1][start_col-1] + algoScores.mismatchScore;
	n2 = scoring_matrix[start_row][start_col-1] + algoScores.gapScore;  //left
	n3 = scoring_matrix[start_row-1][start_col] + algoScores.gapScore;  //up

	//check if this cell have same value with diagonally left cell
	if(scoring_matrix[start_row][start_col]==n1){  //diagonally
		results_vec[index].result_pair.q = sequences.q[start_row-1] + results_vec[index].result_pair.q;
		results_vec[index].result_pair.d = sequences.d[start_col-1] + results_vec[index].result_pair.d;
		//if the diagonally left cell is zero then we are at first row,column,so end of recursion
		if(scoring_matrix[start_row-1][start_col-1]==0){
			results_vec[index].start = start_col-1;
			return steps;
		}
		else//continue recursion
			return traceback(steps+1, start_row-1, start_col-1, sequences, scoring_matrix, results_vec);
	}
	//check if this cell have same value with left cell
	else if(scoring_matrix[start_row][start_col]==n2){ //left
		results_vec[index].result_pair.q = "-" + results_vec[index].result_pair.q;
		results_vec[index].result_pair.d = sequences.d[start_col-1] + results_vec[index].result_pair.d;
		//if the  left cell is zero then we are at first column, so end of recursion
		if(scoring_matrix[start_row][start_col-1]==0){
			results_vec[index].start = start_col-1;
			return steps;
		}
		else//continue recursion
			return traceback(steps+1, start_row, start_col-1, sequences, scoring_matrix, results_vec);
	}//check if this cell have same value with up cell
	else if(scoring_matrix[start_row][start_col]==n3){ //up
		results_vec[index].result_pair.q = sequences.q[start_row-1] + results_vec[index].result_pair.q;
		results_vec[index].result_pair.d = "-" + results_vec[index].result_pair.d;
		//if the  up cell is zero then we are at first row,so end of recursion
		if(scoring_matrix[start_row-1][start_col]==0){
			results_vec[index].start = start_col-1;
			return steps;
		}
		else//continue recursion
			return traceback(steps+1, start_row-1, start_col, sequences, scoring_matrix, results_vec);
	}
	return 0;
}
