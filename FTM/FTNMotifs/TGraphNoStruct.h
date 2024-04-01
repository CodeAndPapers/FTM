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
		delete[]EMaxIntvl;
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
		int choiceStartT, int choiceEndT);

	void findTMotifsDynamic(int k, vec(TMotif*)*& newResult, int oriEndT,
		i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber) {}

//	void genMotifInOneIntv(SAVEINFONOST_VIter& iterStart, SAVEINFONOST_VIter& iterEnd,
//		i2iHMap& vertex2Pos, DisjointSet*& disjointSet, int& vertexNum,
//		vec(int)& combineCCPos, int& realMotifNum,
//		i2iHMap& root2Comp, vec(CComponents*)& tempComponents,
//		vec(int)& saveCCPos, int motifStartT, int motifEndT,
//		vec(TMotif*)*& result, int k, long long& motifNumber);
//
//
//#pragma region generateMaxTM
//	/*fetch the new edge from one R edge set and insert into the disjoint set (maintain connectivity)*/
//	void maintainUFSet(SAVEINFONOST_VIter& infoBegin,
//		SAVEINFONOST_VIter& infoEnd, i2iHMap& vertex2Pos, DisjointSet*& ufset, int&vertexNum,
//		vec(int)& combineCCPos);
//
//	/*add edge into generated motif or generate new motif and check left expandable
//			(generateMaxTM case 1, 2 and add edge into generated motif in generateMaxTM case 3)*/
//	void updateNewEdgeInfo(
//		SAVEINFONOST_VIter& infoBegin, SAVEINFONOST_VIter& infoEnd,
//		vec(CComponents*)& tempComponents,
//		i2iHMap& vertex2Pos, DisjointSet*& disjointSet,
//		i2iHMap& root2Comp, /*int& tempComponentsSize,*/ int& realMotifNum,
//		vec(int)& saveCCPos, int startTime);
//#pragma endregion 
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

		
	void computeRES(int intvB, int intvE,
		vec(int)*& edgeSetsR, int& selectedNum,
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed) {}
	/*computeRES for FTM*/
	void computeRES(int choiceStartT, int choiceEndT,
		SAVEINFONOST_Vec*& edgeSetsR,
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed, int k, int setsRNum);
	void computeRES(int intvB, int intvE, int*& scanP,
		SAVEINFO_Vec*& edgeSetsR, int& selectedNum,
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed);

	/*computeRES for DFTM (row <= T-k+1)*/
	void computeRESForDFTM(int intvB, int intvE,
		vec(int)*& edgeSetsR, int& selectedNum,
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed, int k, vec(TMotif*)*& result, int pos) {}

	/*get the label of edge with edgeId at the time*/
	Label getWeight(int time, int edgeId);

#pragma region structure 
#pragma region Only save temporal graph label   
	Label** eW; // lab_t:temporal graph label   eW[t][edgeId] O(TE)
#pragma endregion
#pragma endregion  
};