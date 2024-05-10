#pragma once
#include"stdafx.h"
#include"TGraph.h"
#include"ICTree.h"


/*use IC trees*/
class TGraphUICTree :public TGraph {
public:

	TGraphUICTree() { cout << "choose new tree" << endl; }
	// copy-constructed
	TGraphUICTree(const TGraphUICTree& ances);

	//output the information of TGraph conveniently
	friend ostream& operator<<(ostream& output, const TGraphUICTree& e) {
		output << "Input temporal graph information:" << endl;
		output << "number of node: " << e.nNode << "\tnumber of edge: " << e.nEdge << endl;
		output << "time length: " << e.nTimestamp <<
			"\tstart time: " << e.startT << "\tend time: " << e.endT << endl;
		return output;
	}

	~TGraphUICTree() {
		for (int i = 0; i < nEdge; i++) {
			delete tree[i];
		}
		delete[]tree;
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
	void changeGraph(const char* src, int endT, int limitNewEndT);
	#pragma endregion

	void findTMotifs(int k, vec(TMotif*)*& result,
		i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber,
		int choiceStartT, int choiceEndT, int TFchoice);
private:
	/*
	create the eL table and IC trees
	four vectors save four number u,v,t,w respectively
		u,v:two endpoints of an edge
		t: timestamp
		w: edge label
	*/
	void createStructure(vec(int)& u_arr, vec(int)& v_arr,
		vec(int)& t_arr, vec(Label)& w_arr);
	
	/*computeRES for FTM*/
	void computeRES(int intvB, int intvE,
		vec(int)*& edgeSetsR, int& selectedNum,
		unordered_map<int, bool>& fixLabel,
		bool isEdgeTypeFixed);
	
	/*computeRES for DFTM (row number<=T-k+1)*/
	void computeRESForDFTM(int intvB, int intvE,
		vec(int)*& edgeSetsR, int& selectedNum,
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed, int k, vec(TMotif*)*& result, int pos);

	/*get the label of edge with edgeId at the time*/
	Label getWeight(int time, int edgeId);

	#pragma region auxiliary structure 
	ICTree** tree;// trees for every edge
	#pragma endregion
	
};