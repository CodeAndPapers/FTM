#pragma once
#include "stdafx.h"
#include "TGraph.h"
#include "TMotif.h"
#include "Interval.h"
#include "Util.h"
#include "DisjointSet.h"

/*implement algorithms*/
class FindTMotif {
public:
	using Intv = Interval<int, Label>;

	#pragma region FTM
	/*
		@parameter:
		graph: input temporal graph
		result: output result TF
	*/
	static void FTM(TGraph*& graph,
		 vec(TMotif*)*& result, 
		i2bHMap& fixLabel/*, int methodId*/);
	#pragma endregion 

	#pragma region DFTM
	/*
		@parameter:
		graph: input temporal graph (updated)
		result: original result TF^o
		newResult: new result TF
		oriEndT: original ending timestamp of temporal graph
	*/
	static void DFTM(TGraph*& graph,
		vec(TMotif*)*& result, vec(TMotif*)*& newResult,
		int oriEndT, i2bHMap& fixLabel/*, int methodId*/);
	#pragma endregion
	
	//inialization
	/*static void preprocess(TGraph*& graph) {
		graph->resetStructure();
	}*/
	static bool cmp(TMotif*& a, TMotif*& b) {
		sort(a->getMotifEdge()->begin(), a->getMotifEdge()->end());
		sort(b->getMotifEdge()->begin(), b->getMotifEdge()->end());
		return (a->getSize()<b->getSize())||
			(a->getSize() == b->getSize()&&a->getMotifEdge()->begin()->id < b->getMotifEdge()->begin()->id);
	}
	/*print the list of motifs*/
	static void print(TGraph*& temporal_graph, vec(TMotif*)*& res,int pos, bool outputTime) {
		vec(TMotif*) lis = res[pos];
		if (lis.size() == 0) {
			cout << "empty result\n" << endl;
			return;
		}
		sort(lis.begin(), lis.end(), cmp);
		veciter(TMotif*) listIter = lis.begin(), 
			listEnd = lis.end();
		if (listIter != listEnd) {
			int intvLen = (*listIter)->getEndT()- (*listIter)->getStartT() + 1;
			Test::maxIntvLen = Test::maxIntvLen < intvLen ? intvLen : Test::maxIntvLen;
			Test::sumIntvLen += lis.size() * intvLen;
		}
		int i = 1;
		TMotif* motif;
		int size;
		for (; listIter != listEnd; ++listIter) {
			motif = *listIter;

			size = motif->getSize();
			

			motifSum += size;
			motifMaxNum = size <= motifMaxNum ? motifMaxNum : size;
			motifMinNum = size >= motifMinNum ? motifMinNum : size;
			/*if (output == 3) {
				if (Test::hist.find(size) == Test::hist.end()) Test::hist[size] = 1;
					else Test::hist[size]++;
				temporal_graph->checkVertex(motif); 
			}*/
			if (output == 1||output == 2) {
				if (outputTime) {
					cout << "startT: " << motif->getStartT()
						<< "\tendT: " << motif->getEndT() << endl;
				}
				if (output == 2) {
					cout << MOTIF_ID << i++ << endl;
					cout << EDGE_NUM << size << endl;
					if (size > 100) continue;
					temporal_graph->printMotif(res, motif, k);
				}
				else if (output == 1) {
					cout << MOTIF_ID << i++ << endl;
					cout << EDGE_NUM << size << endl;
				}
			}
		}
		if (output != 0) cout << "\n";
	}

	static int k;//frequency condition

	static long long motifNumber;//the number of motifs

	static int output;//output level
						   
	static bool isEdgeTypeFixed;//whether fix the edge's label of motifs(default:false)

	static char outputSrc[FILE_NAME_LENGTH];//output file name

	static long long  motifMaxNum, motifMinNum, motifSum;
};