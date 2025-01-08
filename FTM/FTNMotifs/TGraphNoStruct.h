#pragma once
#include "stdafx.h"
#include "TGraph.h"

/*use EL table*/
class TGraphNoStruct :public TGraph {
public:

	TGraphNoStruct() { cout << "choose no create struct" << endl; }
	// copy-constructed
	TGraphNoStruct(const TGraphNoStruct& ances);

	//output the information of TGraph conveniently
	friend ostream& operator<<(ostream& output, const TGraphNoStruct& e) {
		output << "Input temporal graph information:" << endl;
		output << "number of node: " << e.nNode << "\tnumber of edge: " << e.nEdge << endl;
		output << "time length: " << e.nTimestamp <<
			"\tstart time: " << e.startT << "\tend time: " << e.endT << endl;
		return output;
	}

	~TGraphNoStruct() {
		for (int i = 0; i < nTimestamp; i++) {
			delete[]eW[i];
		}
		delete[]eW; 
		delete[]cE;
	}

#pragma region construct and update the temporal graph
	/*
	construct the temporal graph
	1) graph file: each line of the file describe an edge in a timestamp.
	each line has four number: u,v,t,w (separated by ',' with no other space)
	for weight w of edge (u,v) in time t. Ids of node u and v are not guaranteed to be
	continuous, while the timestamps t are continuous, i.e. all edges in time 0 come
	first, and followed by edges in time 1, 2...
	*/
	void constructGraph(const char* src);

	/*construct the temporal graph, but fixed the startT and endT*/
	void constructGraph(const char* src, int startT, int endT);

	//increase snapshots of the temporal graph
	void changeGraph(const char* src, int oldEndT, int limitNewEndT);
#pragma endregion

	void findTMotifs(int k, vec(TMotif*)*& result,
		i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber,
		int choiceStartT, int choiceEndT, int TFchoice);
	void findTMotifsDynamic(int k, vec(TMotif*)*& newResult, int oriEndT,
		i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber) {}

	/*add edge into generated motif or generate new motif and check left expandable
			(generateMaxTM case 1, 2 and add edge into generated motif in generateMaxTM case 3)*/
	void TGraphNoStruct::updateNewEdgeInfo(
		veciter(int)& infoBegin, veciter(int)& infoEnd,
		vec(CComponents*)& tempComponents,
		i2iHMap& vertex2Pos, DisjointSet*& disjointSet,
		i2iHMap& root2Comp, /*int& tempComponentsSize,*/ int& realMotifNum,
		vec(int)& saveCCPos, int startTime);
	/*DFTM (row number<=T-k+1)*/
	void findTMotifsDynamic(int k, vec(TMotif*)*& newResult, int oriEndT,
		i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber, int TFchoice);

	void runDFTM(int k, vec(TMotif*)*& newResult, int oriEndT, i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber, int TFchoice);
private:
	/*
	create the eL table and IC trees
	four vectors save four number u,v,t,w respectively
		u,v:two endpoints of an edge
		t: timestamp
		w: edge label
	*/
	void createStructure(vec(int)& u_arr, vec(int)& v_arr,
		vec(int)& t_arr, vec(Label)& w_arr/*,
		unordered_map<int, int>* &name2id*/);

	void computeRES(int scanT, int k,
		vec(int)*& edgeSetsR, int& selectedNum,
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed);

	/*computeRES for DFTM (row <= T-k+1)*/
	void computeRESForDFTM(int scanT, int oriEndT, int k,
		vec(int)*& edgeSetsR, int& selectedNum,
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed, vec(TMotif*)*& result, int pos);

	/*get the label of edge with edgeId at the time*/
	Label getWeight(int time, int edgeId);

#pragma region structure 
#pragma region Only save temporal graph label   
	int *cE; // record the ending timestamp of the edge with the same label 
	Label** eW; // lab_t:temporal graph label   eW[t][edgeId] O(TE)
#pragma endregion
#pragma endregion  
};