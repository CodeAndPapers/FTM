#include"stdafx.h"
#include"Util.h"
#include"TGraph.h"
#include"TGraphUEL.h"
#include"TGraphUICTree.h"
#include"TGraphNoStruct.h"
#include"FindTMotif.h"
#include"TMotif.h"

#define LESSMEMFOREUR 1 //method to reduce memory for EURData

//unordered_map<int, int> Test::hist, Test::vertexHist;
clock_t Test::compr, Test::gm, Test::gne;
SIZE_T Test::peakMemory = 0;
long long FindTMotif::motifMaxNum = 0;
long long FindTMotif::motifMinNum = 0x7fffffff;
long long FindTMotif::motifSum = 0; 
long long Test::sumIntvLen, Test::maxIntvLen;
#pragma region function declaration
#pragma region load temporal graph
void createTGraph(TGraph*& temporal_graph,const char* inputSrc, int indexId, int startT, int endT);
#pragma endregion

#pragma region which problem to run
void runStaticAlgorithm(TGraph*& temporal_graph, unordered_map<int, bool>& fixLabel/*,int choice*/ );
void runIncrementalAlgorithm(TGraph*& temporal_graph, const char * src
	, int graphEndT, unordered_map<int, bool>& fixLabel, int newEndT/*, int choice*/);
#pragma endregion

#pragma region read k
void readK(vec(int)& arr, const char* file);
void readFixedLabel(unordered_map<int, bool>& fixLabel,
	const char* file);
#pragma endregion

#pragma region release memory
void releaseResult(vec(TMotif*)*& result, int len);
void releaseResultNotDelete(vec(TMotif*)*& result, int len);
#pragma endregion
#pragma endregion 
/////////////////////////////////////////////////////////////////////////////
/*Input Parameters:
- i : File name of input temporal graph(default: example-temporal-graph.txt)
- f : File name of output(default: example-output.txt)
- k : File name of frequency threshold(default:example-k.txt,if not exists k=10)
		the program can test different frequent conditions for one time(untested)
	for static algorithm:
		file content format:
			10
			20
			30
- r : Algorithm Id (default:1)
		4:static algorithm FTM
		9:dynamic algorithm DFTM
- o : Output level of motif (default:0)
		three levels:
			0:only output motif number, the running time and memory use
			1:except for those outputs mentioned above,
				output the edges number of every motif
			2:except for those outputs mentioned above,
				output the detailed information of motif edges(nodes and labels)
- e : File name of fixed labels of motif edges (default: do not fix labels)
		file content format:
			1
			2 
			(means motifs with label types set = {1,2}, all edges have label type in {1,2},
			one line for one label type)
- g : Index id (default:1)
        1:eL table
		2:interval tree
		5:IC tree
		6:None extra struct
- l : Limit the start time and end time of input data
		(default:doesn't limit) format:-l:0,500
		e.g.: if the interval of input temporal graph is [1,100]
		-l:1,50  means you only use the interval [1,50] 
			of temporal graph when testing the algorithm  (used in dynamic algorithm)
- n : Limit the end time of input data when snapshot increasing
		(use - l at the same time) (default:doesn't limit) 
		format:-l:0,500 -n:1600 
		 means that the interval of temporal graph is [0,500] before snapshots increase,
		 and that the interval of temporal graph is [0,1600] after snapshots increase (used in dynamic algorithm)
- m : Method to save the set TF (default:2)
        1 : maintain the temporal motifs in intervals (the output are in the time order in one row of TI-Table; need O(T^2) space; we use this method for better outputs)
		2 : maintain the temporal motifs in rows (the output are not in the time order in one row of TI-Table; does not need O(T^2) space) 
e.g.:
	for static algorithm FTM:
		(use EL table, only output the number of motifs and edges, time, memory)
		-i:temporalgraph.txt -k:k10.txt -r:4 -g:1 -f:output.txt -o:1
		(use IC trees, output nodes of edges, labels, the number of motifs and edges, time, memory)
		-i:temporalgraph.txt -k:k10.txt -r:4 -g:5 -f:output.txt -o:2
	for dynamic algorithm DFTM:
		(use EL table, only output the number of motifs, time, memory, interval of temporal graph is from [0,500] to [0,1000])
		-i:temporalgraph.txt -k:k10.txt -r:9 -g:1 -f:output.txt -o:0 -l:0,500 -n:1000
*/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
	
	#pragma region parameter settings
	char inputSrc[FILE_NAME_LENGTH] = INPUT_FILE;//input temporal graph file
	char kFile[FILE_NAME_LENGTH] = K_FILE;//file of frequent conditon k
	FindTMotif::k = DEFAULT_K;//frequent condition
	int choice = ALGORITHM_ID;//choose which algorithm to run
	FindTMotif::TFchoice = 2;//choose which method to save the set TF
	char fixedLabelFile[FILE_NAME_LENGTH]
		= FIXEDLABEL_FILE;//file of fixed labels
	int indexId = INDEXID;//index id
	//int runTimes = RUNTIMES;//the testing times
	int startT = -1, endT = -1;//start time, end time of the temporal graph
	int limitEndT = -1;//limit end time when running special case 2
	int newEndT = -1;// new end time of the temporal graph(after snapshot increasing)
	vec(int) fixedK;// multiple frequent conditions
	int sep;//separation of input file 
	FindTMotif::isEdgeTypeFixed = false;//whether the edge label is fixed
	unordered_map<int, bool> fixedLabel;//record fixed labels for one edge
	#pragma endregion

	#pragma region load parameters
		for (int i = 1; i < argc; i++) {
			switch (argv[i][1]) {
				case 'e'://Fixed the label of edge
					strcpy(fixedLabelFile, argv[i] + 3);
					FindTMotif::isEdgeTypeFixed = true;
					break;
				case 'f'://output file name
					strcpy(FindTMotif::outputSrc, argv[i] + 3);
					break;
				case 'g'://choose index in TGraph
					indexId = STR2INT(argv[i] + 3);
					break;
				case 'i'://File name of input graph
					strcpy(inputSrc,argv[i] + 3);
					break;
				case 'k'://Frequency threshold
					strcpy(kFile, argv[i] + 3);
					break;
				case 'l'://fixed startT and endT
					sep = (int)(find(argv[i] + 3, argv[i] + 3 + LINE_LENGTH, SEP_CHAR) - (argv[i] + 3));
					startT = STR2INT(argv[i] + 3);
					endT = STR2INT(argv[i] + 4 + sep);
					break;
				case 'n'://fixed newEndT
					newEndT = STR2INT(argv[i] + 3);
					break;
				case 'o'://output the motif or not
					FindTMotif::output = STR2INT(argv[i] + 3);
					break;
				case 'r'://Choose Algorithm Id
					choice = STR2INT(argv[i] + 3);
					break;
				case 'm'://Choose saving method for the set TF
					FindTMotif::TFchoice = STR2INT(argv[i] + 3);
					break;
				default://do nothing
					break;
			}
		}
	
		//argument setting
		OUTPUT_FILE_OPEN
		cout << "graph : " << inputSrc << endl;
		readK(fixedK, kFile);
		cout << "choice : " << choice << endl;
		if (startT != -1) {
			cout << "fixed time : " << startT << "," << endT << endl;
		}
	#pragma endregion

	#pragma region load fixed label
		if (FindTMotif::isEdgeTypeFixed) {
			readFixedLabel(fixedLabel, fixedLabelFile);
		}
	#pragma endregion

	#pragma region load graph
	TGraph* temporal_graph;
	createTGraph(temporal_graph, inputSrc, indexId, startT, endT);
	#pragma endregion

	#pragma region run the algorithm
		switch (choice) {
			case 4://static algorithm FTM
				for (auto iter = fixedK.begin(); iter != fixedK.end(); ++iter) {
					FindTMotif::k = *iter;
					cout << "k : " << *iter << endl;
					runStaticAlgorithm(temporal_graph,
						fixedLabel/*, choice*/);
				}
				FindTMotif::output = 0;
				break;
			case 9: //dynamic algorithm DFTM
				for (auto iter = fixedK.begin(); iter != fixedK.end();
						++iter) {
					FindTMotif::k = *iter;
					cout << "k : " << *iter << endl;
					runIncrementalAlgorithm(
						temporal_graph, inputSrc,
						endT, fixedLabel, newEndT/*, choice*/);
				}
				FindTMotif::output = 0;
				break;
			default:
				cout << "please choose a algorithm to run" << endl;
				return 0;
		}
	#pragma endregion
	
	delete temporal_graph;
	cout << "\n";

	_CrtDumpMemoryLeaks();//check memory leak
	OUTPUT_FILE_CLOSE
}

#pragma region function implementation
#pragma region load temporal graph
void createTGraph(TGraph*& temporal_graph,const char* inputSrc, int indexId,
		int startT, int endT) {
	switch (indexId) {
		case 1://eL table
			temporal_graph = DBG_NEW TGraphUEL();
			break;
		//case 2://interval tree
		//	temporal_graph = DBG_NEW TGraphUIntvTree();
		//	break;
		case 5://IC trees
			temporal_graph = DBG_NEW TGraphUICTree();
			break;
		case 6://no struct
			temporal_graph = DBG_NEW TGraphNoStruct();
			break;
		default:
			cout << "wrong index id" << endl;
			exit(-1);
	}
	if (startT != -1) {
		temporal_graph->constructGraph(inputSrc, startT, endT);
	}
	else {
		temporal_graph->constructGraph(inputSrc);
	}
	
	cout << *temporal_graph; 
	cout << "after create index: "; Test::showMemoryUse();
	//if (indexId == 6) exit(0);
}
#pragma endregion

#pragma region which problem to run

#pragma region static algorithm
void runStaticAlgorithm(TGraph*& temporal_graph, unordered_map<int, bool>& fixLabel/*, int choice*/) {
	#pragma region initialzation
	int resultLen;
	vec(TMotif*)* lis=NULL;
	size_t size;
	FindTMotif::motifNumber = 0;
	int nTimestamps = temporal_graph->getNTimestamp();
	if (FindTMotif::TFchoice == 2) {
		resultLen = (nTimestamps - FindTMotif::k + 1) << 1;
	}
	else {
		resultLen = (nTimestamps - FindTMotif::k + 2)*(nTimestamps - FindTMotif::k + 1) >> 1;
	}
	lis = DBG_NEW vec(TMotif*)[resultLen];//final result
	#pragma endregion
	clock_t startTime, endTime;
	Test::gne = Test::gm = Test::compr = 0;
	Test::peakMemory = 0;
	cout << "before algorithm: "; Test::showMemoryUse();
	startTime = clock();
	FindTMotif::FTM(temporal_graph, lis, fixLabel/*,
		choice==3 ? 1 : 2*/);
	endTime = clock();
	#pragma region output
	//cout << "Max Cost:" << Test::maxCost << endl;
	//cout << "All Cost(before optimized):" << Test::allCost1 << endl;
	//cout << "All Cost(after optimized):" << Test::allCost2 << endl;
	//cout << "edges number:" << Test::edgesNum << endl;
	cout<<"compr: "<< Test::compr << "ms, gm: " << Test::gm << "ms, gne: "<< Test::gne << "ms" << endl;
	cout << "time: " << endTime - startTime << "ms" << endl;
	veciter(TMotif*) iter;
	TMotif* motif;
	FindTMotif::motifMaxNum = FindTMotif::motifSum = 0;
	FindTMotif::motifMinNum = 0x7fffffff;
	Test::maxIntvLen = 0, Test::sumIntvLen = 0;
	if (FindTMotif::TFchoice == 2) {
		for (int i = 0; i < resultLen; i += 2) {
			size = lis[i].size() + lis[i + 1].size();//one row 
			if (size == 0) continue;

			if (FindTMotif::output >= 1) {
				cout << MOTIF_NUM << size << endl;
				iter = lis[i].begin();
				motif = *iter;
				//cout << "startT: " << motif->getStartT()
				//	<< "\tendT: " << motif->getEndT() << endl;
			}

			FindTMotif::print2(temporal_graph,
				lis, i, false);
		}
	}
	else {
		for (int i = 0; i < resultLen; i++) {
			size = lis[i].size();
			if (size == 0) continue;

			if (FindTMotif::output >= 1) {
				cout << MOTIF_NUM << size << endl;
				iter = lis[i].begin();
				motif = *iter;
				cout << "startT: " << motif->getStartT()
					<< "\tendT: " << motif->getEndT() << endl;
			}

			FindTMotif::print(temporal_graph,
				lis, i, false);
		}
	}
	
	cout << "motif max edges num: " << FindTMotif::motifMaxNum << " motif min edges num: " << FindTMotif::motifMinNum <<
		"\tmotif avg edges num: " << FindTMotif::motifSum * 1.0 / FindTMotif::motifNumber << endl;
	cout << "motif max intverval length: " << Test::maxIntvLen <<
		"\tmotif avg intverval length: " << Test::sumIntvLen * 1.0 / FindTMotif::motifNumber << endl;
	cout << "sum: " << FindTMotif::motifNumber << endl;
	cout << "after algorithm: "; Test::showPeakMemoryUse();
	Test::showRealPeakMemoryUse();
	//Test::showMemoryUse();
	
	#pragma endregion
	releaseResult(lis, resultLen);//release the memory
	cout << "release result: ";
	Test::showMemoryUse(); cout << "\n";
}
#pragma endregion

#pragma region dynamic algorithm
void runIncrementalAlgorithm(TGraph*& temporal_graph, const char * src
	, int graphEndT, unordered_map<int, bool>& fixLabel, int newEndT/*, int choice*/) {
	if (graphEndT == -1)graphEndT = temporal_graph->getEndT();
	clock_t startTime, endTime;
	int nTimestamp = temporal_graph->getNTimestamp();
	int listSize = nTimestamp - FindTMotif::k + 1;
	int resultLen;
	if (FindTMotif::TFchoice == 2) {
		resultLen = listSize << 1;//original result size
	}
	else {
		resultLen = (listSize + 1)*(listSize) >> 1;//original result size
	}
	int newResultLen;//new result size
	vec(TMotif*)* result = DBG_NEW vec(TMotif*)[resultLen], //original result of general case
		*newResult = NULL;//new result of general case
		//,*midResult = NULL//intermediate result of general case
	size_t size; 
	#pragma region get original result of general case
	startTime = clock();
	FindTMotif::FTM(temporal_graph, result, fixLabel/*, 2*/);
	endTime = clock();
	cout << "gctime: " << endTime - startTime << "ms" << endl;
	#pragma endregion

	#pragma region get Etf result size
	int startT = temporal_graph->getStartT(), endT = graphEndT;
	int allEdgeSize = 0, maxEdgeSize = 0;
	int lastRow = endT - FindTMotif::k + 1;
	for (int i = startT; i <= lastRow; i++) {
		int tempPos;
		if (FindTMotif::TFchoice == 2) {
			tempPos = ((i - startT) << 1) + 1;
		}
		else {
			tempPos = resultPos(i, graphEndT, startT, graphEndT, FindTMotif::k);
		}
		auto resultEnd = result[tempPos].end();
		size = 0;
		for (auto resultIter = result[tempPos].begin(); resultIter != resultEnd; ++resultIter) {
			TMotif* motif = (TMotif*) *resultIter;
			size += motif->getSize();
		}
		allEdgeSize += size;
		maxEdgeSize = maxEdgeSize < size ? size : maxEdgeSize;
	}
	#pragma endregion
	cout << "averEtf: "  << allEdgeSize/(lastRow - startT + 1) << " maxEtf: " << maxEdgeSize << endl;
	#pragma region snapshots increase
	//delete temporal_graph;
	//createTGraph(temporal_graph, src, indexofId, startT, newEndT);
	temporal_graph->changeGraph(src, graphEndT, newEndT);
	cout << *temporal_graph;
	#pragma endregion
	#pragma region preprocess
	FindTMotif::motifNumber = 0;
	listSize = temporal_graph->getNTimestamp() - FindTMotif::k + 1;
	endT = temporal_graph->getEndT();
	/*size of new intermediate result after snapshots increasing*/
	if (FindTMotif::TFchoice == 2) {
		newResultLen = listSize << 1;
	}
	else {
		newResultLen = (listSize + 1)*listSize >> 1;
	}
	newResult = DBG_NEW vec(TMotif*)[newResultLen];//new result
	Test::peakMemory = 0;
	cout << "before algorithm: "; Test::showMemoryUse();
	#pragma endregion
	//int methodid = choice == 8 ? 1 : 2;
	startTime = clock();
	FindTMotif::DFTM
	(temporal_graph, 
		result, newResult, graphEndT, fixLabel/*, methodid*/);
	endTime = clock();
	#pragma region output
	cout << "inctime: " << endTime - startTime << "ms" << endl;
	veciter(TMotif*) iter;
	TMotif* motif;
	if (FindTMotif::TFchoice == 2) {
		for (int i = 0; i < newResultLen; i += 2) {
			size = newResult[i].size() + newResult[i + 1].size();
			if (size == 0) continue;
			if (FindTMotif::output >= 1) {
				cout << MOTIF_NUM << size << endl;
				//iter = newResult[i].begin();
				//motif = *iter;
				//cout << "startT: " << motif->getStartT()
				//	<< "\tendT: " << motif->getEndT() << endl;
			}

			FindTMotif::print2(temporal_graph,
				newResult, i, false);
		}
	}
	else {
		for (int i = 0; i < newResultLen; i++) {
			size = newResult[i].size();
			if (size == 0) continue;
			if (FindTMotif::output >= 1) {
				cout << MOTIF_NUM << size << endl;
				iter = newResult[i].begin();
				motif = *iter;
				cout << "startT: " << motif->getStartT()
					<< "\tendT: " << motif->getEndT() << endl;
			}

			FindTMotif::print(temporal_graph,
				newResult, i, false);

		}
	}
	cout << "motif max edges num: " << FindTMotif::motifMaxNum <<
		"\tmotif avg edges num: " << FindTMotif::motifSum * 1.0 / FindTMotif::motifNumber << endl;
	cout << "sum: " << FindTMotif::motifNumber << endl;
	cout << "after algorithm: "; Test::showPeakMemoryUse();
	Test::showRealPeakMemoryUse();
	//Test::showMemoryUse();
	#pragma endregion 
	
	//releaseResultNotDelete(result, resultLen);
	releaseResult(newResult, newResultLen);
	
	cout << "release result: ";
	Test::showMemoryUse(); cout << "\n";
}
#pragma endregion
#pragma endregion 

#pragma region read k and fixed label from files
//load the setting of frequent condition k
void readK(vec(int)& arr, const char* file) {
	FILE* f;
	f = fopen(file, "r+");
	if (!f) {
		arr.emplace_back(DEFAULT_K);
	}
	else {
		char line[LINE_LENGTH];
		CLEARALL(line, 0, LINE_LENGTH, char);
		while (fgets(line, LINE_LENGTH, f)) {
			//FindTMotif::motifNumber = 0;
			if (strlen(line) == 0) continue;
			arr.emplace_back(STR2INT(line));
		}
		fclose(f);
	}
}

 /*load the setting of fixed labels */
void readFixedLabel(unordered_map<int, bool>& fixLabel,
	const char* file) {
	FILE* f;
	f = fopen(file, "r+");
	if (!f) {
		FindTMotif::isEdgeTypeFixed = false;
		return;
	}
	else {
		char line[LINE_LENGTH];
		CLEARALL(line, 0, LINE_LENGTH, char);
		while (fgets(line, LINE_LENGTH, f)) {
			if (strlen(line) == 0) continue;
			fixLabel[STR2INT(line)] = true;
		}
		fclose(f);
	}
}

#pragma endregion

#pragma region release memory
void releaseResult(vec(TMotif*)*& result, int len) {
	size_t size;
	for (int i = 0; i < len; i++) {
		size = result[i].size();
		for (size_t j = 0; j < size; j++) {
			delete result[i][j];
		}
		result[i].clear();
	}
	delete[] result;
}
void releaseResultNotDelete(vec(TMotif*)*& result, int len) {
	for (int i = 0; i < len; i++) {
		result[i].clear();
	}
	delete[] result;
}
#pragma endregion
#pragma endregion 