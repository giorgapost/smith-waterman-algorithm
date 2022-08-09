/**
 * @file Framework.cpp
 */
#include "Framework.h"

Framework::Framework(int argc, char* argv[]) {
    if(argc<=1)
        readArgsFromConsole();
    else{ //the user provided the arguments, so initialize with them
        bool argName = false, argPath = false, argMatch = false, argMismatch = false, argGap = false;
        for(int i=0;i<argc-1;i+=1){
            if(!string(argv[i]).compare("-id")){
                reportId = string(argv[i+1]);
                argName = true;
            }
            else if(!string(argv[i]).compare("-path")){
                path = string(argv[i+1]);
                argPath = true;
            }
            else if(!string(argv[i]).compare("-match")){
                algoScores.matchScore = atoi(argv[i+1]);
                argMatch = true;
            }
            else if(!string(argv[i]).compare("-mismatch")){
                algoScores.mismatchScore = atoi(argv[i+1]);
                argMismatch = true;
            }
            else if(!string(argv[i]).compare("-gap")){
                algoScores.gapScore = atoi(argv[i+1]);
                argGap = true;
            }
            else
                continue;
        }

        if(!argName || !argPath || !argMatch || !argMismatch || !argGap)
            throw invalid_argument("Error. Missing arguments.");
    }

	readInputFile(); //read the sequences for alignment from the input file

	for(size_t i=0;i<data.size();i++) //initialize the results vector
        results.push_back(vector<Result>());

    statisticData.cellsGreaterThanZero = 0; //Initialize the object holding the statistical data
    statisticData.totalTracebackSteps = 0;
    statisticData.totalTime = 0;
    statisticData.calcCellsTime = 0;
    statisticData.totalTracebackTime = 0;
}

Framework::~Framework() {}

long double Framework::getTime(void){
	struct timeval ttime;
	gettimeofday(&ttime, 0);
	return (ttime.tv_sec+ttime.tv_usec*0.000001);
}

void Framework::printResultsToFile(void){
    ofstream output;
    output.open ((string("reports/Report_") + reportId + string(".txt")).c_str());
	if(!output.is_open())
		throw ios_base::failure((string("Output file reports/Report_") + reportId + string(".txt cannot be opened.")).c_str());


    for(size_t i=0;i<results.size();i++){
        output << "Q: " << data[i].q << endl;
        output << "D: " << data[i].d << endl;

        for(size_t j=0;j<results[i].size();j++){
            output << "Match " << j+1 << " [Score: " << results[i][j].score << ", Start: " << results[i][j].start << ", Stop: " << results[i][j].stop << "]" << endl;
            output << "	D: " << results[i][j].result_pair.d << endl;
            output << "	Q: " << results[i][j].result_pair.q << endl;
        }
    }
    output.close();
}

void Framework::readArgsFromConsole(void){
    cin.clear();
    cout << "Please enter the path to the file with the input data:";
    cin >> path;

    cin.clear();
    cout << "Please enter an ID for the report that will be extracted:";
    cin >> reportId;

    cin.clear();
    cout << "Please enter an integer for the match score of the algorithm:";
    cin >> algoScores.matchScore;

    cin.clear();
    cout << "Please enter an integer for the mismatch score of the algorithm:";
    cin >> algoScores.mismatchScore;

    cin.clear();
    cout << "Please enter an integer for the gap score of the algorithm:";
    cin >> algoScores.gapScore;
}

void Framework::readInputFile(void){
	ifstream input;
	string tmp_string;
	Pair qd_pair;

	input.open((this->path).c_str()); //Open the file with the input data
	if(!input.is_open()){
		input.close();
		throw ios_base::failure("Input file cannot be opened.");
	}
	else{  //File opened
		input >> tmp_string;  //"Q:" string of 1st pair
		while(!input.eof()){  //While not reached EOF, continue reading the next pair
			if(tmp_string.compare("Q:")){  //1st sequence of a pair not starting with "Q:"
                input.close();
				throw runtime_error("Error: Invalid file contents...");
			}
			else{  //1st sequence starts with "Q:" - all good
				while(tmp_string.compare("D:")){  //While not found "D:", continue reading
					if(input.eof()){  //EOF found before "D:" --> error
						input.close();
                        throw runtime_error("Error: Invalid file contents...");
					}

					input >> tmp_string;
					if(tmp_string.compare("D:")) //if not found "D:", we are still reading q
						qd_pair.q += tmp_string;
				}
				while(tmp_string.compare("Q:") && !input.eof()){ //While not found EOF or "Q:" of the next pair, continue reading
					tmp_string = "";
					input >> tmp_string;
					if(tmp_string.compare("Q:")) //if not found "Q:", we are still reading d
						qd_pair.d += tmp_string;
				}
				data.push_back(qd_pair); //push the data to the vector
				qd_pair.q.clear();
				qd_pair.d.clear();
			}
		}
		input.close();  //close the ifstream
	}
}
