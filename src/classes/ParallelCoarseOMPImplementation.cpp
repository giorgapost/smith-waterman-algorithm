/**
 * @file ParallelCoarseOMPImplementation.cpp
 */
#include "ParallelCoarseOMPImplementation.h"

ParallelCoarseOMPImplementation::ParallelCoarseOMPImplementation(int argc, char* argv[]): Framework(argc, argv) {
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

    statisticData.calcCellsTime = new long double[threads];
    statisticData.totalTracebackTime = new long double[threads];
    for(size_t i=0;i<threads;i++){
        statisticData.calcCellsTime[i] = 0;
        statisticData.totalTracebackTime[i] = 0;
    }
}

ParallelCoarseOMPImplementation::~ParallelCoarseOMPImplementation() {
    delete statisticData.calcCellsTime;
    delete statisticData.totalTracebackTime;
}

long long int ParallelCoarseOMPImplementation::fill_scoring_matrix(Pair sequences, vector< vector<int> >& scoring_matrix, vector<Position>& max_pos){
	int n1, n2, n3;
	int maxCell=0;
	unsigned int row, col;
	Position tmp_pos;
	long long int cellsGreaterThanZero=0;

	max_pos.clear();

	for(row=0;row<scoring_matrix.size();row++){
		for(col=0;col<scoring_matrix[row].size();col++){
			if((row==0) || (col==0)){
				scoring_matrix[row][col] = 0;
			}
			else{
				if(sequences.q[row-1]==sequences.d[col-1])
					n1 = scoring_matrix[row-1][col-1] + algoScores.matchScore;
				else
					n1 = scoring_matrix[row-1][col-1] + algoScores.mismatchScore;
				n2 = scoring_matrix[row-1][col] + algoScores.gapScore;
				n3 = scoring_matrix[row][col-1] + algoScores.gapScore;
				if(n1>=n2 && n1>=n3){
					if(n1>=0)
						scoring_matrix[row][col] = n1;
					else
						scoring_matrix[row][col] = 0;
				}
				else if(n2>=n1 && n2>=n3){
					if(n2>=0)
						scoring_matrix[row][col] = n2;
					else
						scoring_matrix[row][col] = 0;
				}
				else{  // if(n3>=n1 && n3>=n2)
					if(n3>=0)
						scoring_matrix[row][col] = n3;
					else
						scoring_matrix[row][col] = 0;
				}
			}

			if(scoring_matrix[row][col]>0)
				statisticData.cellsGreaterThanZero++;

			if(scoring_matrix[row][col]>maxCell){
				max_pos.clear();
				tmp_pos.row=row;
				tmp_pos.column=col;
				max_pos.push_back(tmp_pos);
				maxCell = scoring_matrix[row][col];
			}
			else if(scoring_matrix[row][col]==maxCell){
				tmp_pos.row=row;
				tmp_pos.column=col;
				max_pos.push_back(tmp_pos);
			}
		}
	}
	return cellsGreaterThanZero;
}

void ParallelCoarseOMPImplementation::printStatistics(void){
    cout << "A) Total pairs of sequences Q-D: " << data.size() << endl;
    cout << "B) Total cells with value: " << statisticData.cellsGreaterThanZero << endl;
    cout << "C) Total traceback steps: " << statisticData.totalTracebackSteps << endl;
    cout << "D) Total time of program execution: " << statisticData.totalTime << " seconds"<< endl;
    cout << "E) Cell Updates Per Second (CUPS) based on total execution time: " << (double)statisticData.cellsGreaterThanZero/statisticData.totalTime << endl;
    for(size_t i=0;i<threads;i++){
        cout << "F) Thread " << i << " - Total time of calculating cells: " << statisticData.calcCellsTime[i] << " seconds" << endl;
        cout << "G) Thread " << i << " - Total traceback time: " << statisticData.totalTracebackTime[i] << " seconds" << endl;
        cout << "H) Thread " << i << " - Cell Updates Per Second (CUPS) based on time of calculating cells: " << (double)statisticData.cellsGreaterThanZero/statisticData.calcCellsTime[i] << endl;
    }
}

void ParallelCoarseOMPImplementation::runAlgorithm(void){
    long double time0 = getTime();

    #pragma omp parallel shared(data, algoScores, results, statisticData) num_threads(threads)
    {
        Result res;
        vector<Position> max_pos_vec;
        vector<Result> results_vec;
        int cells=0, steps=0;
        double tr_time=0, f_time=0;

        #pragma omp for
        for(size_t i=0;i<data.size();i++){
            vector< vector<int> > scoring_matrix(data[i].q.size()+1, vector<int>(data[i].d.size()+1));

            double time_f1 = getTime();
            cells += fill_scoring_matrix(data[i], scoring_matrix, max_pos_vec);
            double time_f2 = getTime();
            f_time += (time_f2-time_f1);

            for(size_t j=0;j<max_pos_vec.size();j++){
                res.score = scoring_matrix[max_pos_vec[j].row][max_pos_vec[j].column];
                res.stop = max_pos_vec[j].column-1;
                res.ref_id = i;
                results_vec.push_back(res);

                double time_tr1 = getTime();
                steps += traceback(1, max_pos_vec[j].row, max_pos_vec[j].column, data[i], scoring_matrix, results_vec);
                double time_tr2 = getTime();
                tr_time += (time_tr2-time_tr1);
            }
        }
        for(size_t i=0;i<results_vec.size();i++){
            #pragma omp critical
            {
                results[results_vec[i].ref_id].push_back(results_vec[i]);
            }
        }

        #pragma omp critical
        {
            statisticData.calcCellsTime[omp_get_thread_num()] = f_time;
            statisticData.totalTracebackTime[omp_get_thread_num()] = tr_time;
            statisticData.cellsGreaterThanZero += cells;
            statisticData.totalTracebackSteps += steps;
        }
    }
  	statisticData.totalTime = getTime() - time0;
}

long long int ParallelCoarseOMPImplementation::traceback(int steps, int start_row, int start_col, Pair sequences, vector< vector<int> >& scoring_matrix, vector<Result>& results_vec){
	unsigned int index = results_vec.size()-1;
	int n1, n2, n3;


	if(sequences.q[start_row-1]==sequences.d[start_col-1])			 //diagonally up left
		n1 = scoring_matrix[start_row-1][start_col-1] + algoScores.matchScore;
	else
		n1 = scoring_matrix[start_row-1][start_col-1] + algoScores.mismatchScore;

	n2 = scoring_matrix[start_row][start_col-1] + algoScores.gapScore;  //left
	n3 = scoring_matrix[start_row-1][start_col] + algoScores.gapScore;  //up


	if(scoring_matrix[start_row][start_col]==n1){  //diagonally
		results_vec[index].result_pair.q = sequences.q[start_row-1] + results_vec[index].result_pair.q;
		results_vec[index].result_pair.d = sequences.d[start_col-1] + results_vec[index].result_pair.d;

		if(scoring_matrix[start_row-1][start_col-1]==0){
			results_vec[index].start = start_col-1;
			return steps;
		}
		else
			return traceback(steps+1, start_row-1, start_col-1, sequences, scoring_matrix, results_vec);
	}
	else if(scoring_matrix[start_row][start_col]==n2){ //left
		results_vec[index].result_pair.q = "-" + results_vec[index].result_pair.q;
		results_vec[index].result_pair.d = sequences.d[start_col-1] + results_vec[index].result_pair.d;

		if(scoring_matrix[start_row][start_col-1]==0){
			results_vec[index].start = start_col-1;
			return steps;
		}
		else
			return traceback(steps+1, start_row, start_col-1, sequences, scoring_matrix, results_vec);
	}
	else if(scoring_matrix[start_row][start_col]==n3){ //up
		results_vec[index].result_pair.q = sequences.q[start_row-1] + results_vec[index].result_pair.q;
		results_vec[index].result_pair.d = "-" + results_vec[index].result_pair.d;

		if(scoring_matrix[start_row-1][start_col]==0){
			results_vec[index].start = start_col-1;
			return steps;
		}
		else
			return traceback(steps+1, start_row-1, start_col, sequences, scoring_matrix, results_vec);
	}
	return 0;
}
