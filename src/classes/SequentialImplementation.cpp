/**
 * @file SequentialImplementation.cpp
 */
#include "SequentialImplementation.h"

SequentialImplementation::SequentialImplementation(int argc, char* argv[]):
    Framework(argc, argv) {
    statisticData.calcCellsTime = new long double[1];
    statisticData.totalTracebackTime = new long double[1];
}

SequentialImplementation::~SequentialImplementation() {
    delete statisticData.calcCellsTime;
    delete statisticData.totalTracebackTime;
}

long long int SequentialImplementation::fill_scoring_matrix(Pair sequences, vector< vector<int> >& scoring_matrix, vector<Position>& max_pos){
	int n1, n2, n3;
	int max=0;
	Position tmp_pos;
	long long int cellsGreaterThanZero=0;

	max_pos.clear();

	for(size_t row=0;row<scoring_matrix.size();row++){ //for every row of scoring matrix
		for(size_t col=0;col<scoring_matrix[row].size();col++){ //for every column of scoring matrix
			if((row==0) || (col==0)){
				scoring_matrix[row][col] = 0;
			}
			else{
				//compute the cell for match, mismatch and gap
				if(sequences.q[row-1]==sequences.d[col-1])
					n1 = scoring_matrix[row-1][col-1] + algoScores.matchScore;
				else
					n1 = scoring_matrix[row-1][col-1] + algoScores.mismatchScore;
				n2 = scoring_matrix[row-1][col] + algoScores.gapScore;
				n3 = scoring_matrix[row][col-1] + algoScores.gapScore;

				//check which cell has the greater value and fill the cell
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
			//How many cells have value>0
			if(scoring_matrix[row][col]>0)
				cellsGreaterThanZero++;

			//Search for the max value and save the Position x,y of max i matrix
			if(scoring_matrix[row][col]>max){
				max_pos.clear();
				tmp_pos.row=row;
				tmp_pos.column=col;
				max_pos.push_back(tmp_pos);
				max = scoring_matrix[row][col];
			}
			else if(scoring_matrix[row][col]==max){
				tmp_pos.row=row;
				tmp_pos.column=col;
				max_pos.push_back(tmp_pos);
			}
		}
	}
	return cellsGreaterThanZero;
}

void SequentialImplementation::printStatistics(void){
    cout << "A) Total pairs of sequences Q-D: " << data.size() << endl;
    cout << "B) Total cells with value: " << statisticData.cellsGreaterThanZero << endl;
    cout << "C) Total traceback steps: " << statisticData.totalTracebackSteps << endl;
    cout << "D) Total time of program execution: " << statisticData.totalTime << " seconds"<< endl;
    cout << "E) Total time of calculating cells: " << *statisticData.calcCellsTime << " seconds" << endl;
    cout << "F) Total traceback time: " << *statisticData.totalTracebackTime << " seconds" << endl;
    cout << "G) Cell Updates Per Second (CUPS) based on total execution time: " << (double)statisticData.cellsGreaterThanZero/statisticData.totalTime << endl;
    cout << "H) Cell Updates Per Second (CUPS) based on time of calculating cells: " << (double)statisticData.cellsGreaterThanZero/(*statisticData.calcCellsTime) << endl;
}

void SequentialImplementation::runAlgorithm(void){
    Result res;
    long double time0 = getTime();
	for(size_t i=0;i<data.size();i++){
		vector< vector<int> > scoring_matrix(data[i].q.size()+1, vector<int>(data[i].d.size()+1));
        vector<Position> max_pos_vec;
		results[i].clear(); //for the case that this method is accidentally executed more than once

		long double time1 = getTime();
        statisticData.cellsGreaterThanZero += fill_scoring_matrix(data[i], scoring_matrix, max_pos_vec);
		*statisticData.calcCellsTime += getTime() - time1;

		for(size_t j=0;j<max_pos_vec.size();j++){
            res.ref_id = i;
			res.score = scoring_matrix[max_pos_vec[j].row][max_pos_vec[j].column];
			res.stop = max_pos_vec[j].column-1;
			results[i].push_back(res);

            time1 = getTime();
			statisticData.totalTracebackSteps += traceback(1, max_pos_vec[j].row, max_pos_vec[j].column, data[i], scoring_matrix, results[i]);
			*statisticData.totalTracebackTime += getTime() - time1;
		}
	}
	statisticData.totalTime = getTime() - time0;
}

long long int SequentialImplementation::traceback(int steps, int start_row, int start_col, Pair sequences, vector< vector<int> >& scoring_matrix, vector<Result>& results_vec){
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
