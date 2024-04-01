#pragma once
#include "stdafx.h"
#include "TGraph.h"

/*use EL table*/
class TGraphUEL:public TGraph {
public:
	
	TGraphUEL() { cout << "choose eL structure" << endl; }
	// copy-constructed
	TGraphUEL(const TGraphUEL& ances);
	
	//output the information of TGraph conveniently
	friend ostream& operator<<(ostream& output, const TGraphUEL& e) {
		output << "Input temporal graph information:" << endl;
		output << "number of node: " << e.nNode << "\tnumber of edge: " << e.nEdge << endl;
		output << "time length: " << e.nTimestamp <<
			"\tstart time: " << e.startT << "\tend time: " << e.endT << endl;
		return output;
	}

	~TGraphUEL(){
		for (int i = 0; i < nTimestamp; i++) {
			delete[]eL[i];
			delete[]eW[i];
		}
		delete[]eL;
		delete[]eW;
		delete[]maxeL;
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
	void constructGraph(const char* src,int startT,int endT);

	//increase snapshots of the temporal graph
	void changeGraph(const char* src, int oldEndT, int limitNewEndT);
	#pragma endregion

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
	
	/*computeRES for FTM*/
	void computeRES(int intvB, int intvE,
		vec(int)*& edgeSetsR, int& selectedNum,
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed);
	/*void computeRES(int choiceStartT, int choiceEndT,
		SAVEINFO_Vec*& edgeSetsR, 
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed, int k);*/
	
	/*computeRES for DFTM (row <= T-k+1)*/
	void computeRESForDFTM(int intvB, int intvE,
		vec(int)*& edgeSetsR, int& selectedNum,
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed, int k, vec(TMotif*)*& result, int pos);

	/*void computeRESForDFTM(int choiceStartT, int oriEndT, int choiceEndT, int*& scanP,
		SAVEINFO_Vec*& edgeSetsR,
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed, int k, vec(TMotif*)*& result);*/

	/*get the label of edge with edgeId at the time*/
	Label getWeight(int time, int edgeId);

	#pragma region structure 
		#pragma region EL table
			int **eL;// len_t:the times of edges keeping their label fixed O(TE)
			Label** eW; // lab_t:temporal graph label   eW[t][edgeId] O(TE)
		#pragma endregion
		int *maxeL;//max len_t value at t O(T), only used for reducing memory of R edge set
	#pragma endregion  
};