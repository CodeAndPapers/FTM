#include "FindTMotif.h"
#include "TGraphUEL.h"
#include "TGraphUICTree.h"
#include "stdafx.h"

#pragma region parameters initialization
int FindTMotif::k = DEFAULT_K;
int TGraph::maxMotifNum = 0;
long long FindTMotif::motifNumber = 0;
int FindTMotif::output = 0;
int FindTMotif::TFchoice = 0;
bool FindTMotif::isEdgeTypeFixed = false;
char FindTMotif::outputSrc[FILE_NAME_LENGTH] = OUTPUT_FILE;
#pragma endregion 

#pragma region static algorithm
/*  @parameter:
	graph: input temporal graph
	result: output result
*/
void FindTMotif::FTM(TGraph*& graph, vec(TMotif*)*& result,
	i2bHMap& fixLabel/*,int methodId*/ ) {
	int startT = graph->getStartT(), endT = graph->getEndT();
	graph->findTMotifs(k,  result, fixLabel, isEdgeTypeFixed,
		motifNumber, startT , endT, TFchoice);
}
#pragma endregion

#pragma region incremental algorithm
/*
	@parameter:
	graph: input temporal graph (updated)
	result: original result TF^o
	newResult: new result TF
	oriEndT: original ending timestamp of temporal graph
*/
void FindTMotif::DFTM(TGraph*& graph,
	vec(TMotif*)*& result, vec(TMotif*)*& newResult, 
	int oriEndT, i2bHMap& fixLabel/*, int methodId*/) {
	int startT = graph->getStartT();
	int endT = graph->getEndT();
	int nTimestamp = graph->getNTimestamp();
	int nEdge = graph->getNEdge();
	TGraph::maxMotifNum = 0;
	int lastRow = oriEndT - k + 1;
	#pragma region copy from original result
	for (int i = startT; i <= lastRow; i++) {
		int first = i + FindTMotif::k - 1;
		if (TFchoice == 2) {
			int tempPos = (i - startT) << 1, nowPos = (i - startT) << 1;
			int resultSize = (int)result[tempPos].size();
			for (int s = 0; s < resultSize; s++) {
				newResult[nowPos].emplace_back(result[tempPos][s]);
			}
			result[tempPos].clear();
			FindTMotif::motifNumber += resultSize;

			tempPos++;
			nowPos++;//endT
			resultSize = (int)result[tempPos].size();
			for (int s = 0; s < resultSize; s++) {
				newResult[nowPos].emplace_back(result[tempPos][s]);
			}
			if (resultSize > TGraph::maxMotifNum)
				TGraph::maxMotifNum = resultSize;
		}
		else {
			int tempPos = resultPos(i, first, startT, oriEndT, FindTMotif::k) - 1;
			int nowPos = resultPos(i, first, startT, endT, FindTMotif::k) - 1;
			for (int j = first; j <= oriEndT; j++) {
				tempPos++;
				nowPos++;
				int resultSize = (int)result[tempPos].size();
				for (int s = 0; s < resultSize; s++) {
					newResult[nowPos].emplace_back(result[tempPos][s]);
				}
				result[tempPos].clear();
				if (j < oriEndT)
					FindTMotif::motifNumber += resultSize;
				else if (resultSize > TGraph::maxMotifNum)
					TGraph::maxMotifNum = resultSize;
			}
		}
	}
	delete[] result;
	#pragma endregion
	if (oriEndT == endT)return;
	graph->runDFTM(k, newResult, oriEndT, fixLabel, isEdgeTypeFixed, FindTMotif::motifNumber, FindTMotif::TFchoice);
}
#pragma endregion 