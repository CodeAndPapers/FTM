#pragma once
#include "stdafx.h"
#include "TGraph.h"
#include "TMotif.h"
#include "Util.h"
#include "Edge.h"
#include "DisjointSet.h"
#include "ICTree.h"

#define ALLOC_MEM 20000000
/*map the interval [Ts,Te] to the one dimension position of result*/
#define resultPos(Ts,Te,startT,endT,k) ((Ts - startT)*((endT<<1) - startT - Ts - (k<<1) + 5) >>1) + Te - (Ts + k - 1)

/* Temporal Graph: the nodes and edges are fixed
   and the weights vary with time*/
class TGraph {
public:
	using Intv = Interval<int, Label>;

	TGraph() { 
		ufset = NULL;
	}
	// copy-constructed
	TGraph(const TGraph& ances);

	inline int getNTimestamp() const { return nTimestamp; }

	inline size_t getNNode() const { return nNode; }

	inline size_t getNEdge() const { return nEdge; }

	inline int getStartT() const { return startT; }

	inline int getEndT() const { return endT; }

	inline void setStartT(int startT) { this->startT=startT; }

	inline void setEndT(int endT) { this->endT=endT; }

	inline void getEdgeList(Edge*& s) { s = edgeList; }

	inline void getEdge2Ind(map<Edge, int> *& e2i) { e2i = edge2ind; }

	//output the information of TGraph conveniently
	friend ostream& operator<<(ostream& output, const TGraph& e) {
		output << "Input temporal graph information:" << endl;
		output << "number of node: " << e.nNode << "\tnumber of edge: " << e.nEdge << endl;
		output << "time length: " << e.nTimestamp <<
			"\tstart time: " << e.startT << "\tend time: " << e.endT << endl;
		return output;
	}
	
	virtual ~TGraph() {
		delete[] edgeList;
		delete edge2ind;
		if(ufset != NULL)
			delete ufset;
		delete[]EMaxIntvl;
	}

	#pragma region construct and update the temporal graph
	/*
	construct the temporal graph
	graph file src: each line of the file describe an edge in a timestamp.
	each line has four number: u,v,t,w (separated by ',' with no other space)
	for weight w of edge (u,v) in time t. Ids of node u and v are not guaranteed to be
	continuous, while the timestamps t are continuous, i.e. all edges in time 0 come
	first, and followed by edges in time 1, 2...
	*/
	virtual void constructGraph(const char* src)=0;
	
	/*construct the temporal graph, but fixed the startT and endT*/
	virtual void constructGraph(const char* src,int startT,int endT) = 0;

	//increase snapshots of the temporal graph
	virtual void changeGraph(const char* src, int oldEndT, int limitNewEndT) = 0;
	#pragma endregion

	/*print motif*/
	void printMotif(vec(TMotif*)*& res, TMotif* motif, int k, int TFchoice);
	/*count vertex number of motif*/
	//void checkVertex(TMotif* motif);

	#pragma region FTM
		/*FTM*/
	void findTMotifs(int k, vec(TMotif*)*& result,
		i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber,
		int choiceStartT, int choiceEndT, int TFchoice);

	#pragma region FTM step 1
		/*computeRES*/
		virtual void computeRES(int intvB, int intvE,
			vec(int)*& selectedEdge, int& selectedNum,
			unordered_map<Label, bool>& fixLabel,
			bool isEdgeTypeFixed) = 0;
	#pragma endregion 

	#pragma region FTM step 2 and 3
		/*generate motifs in one interval for FTM*/
		void genMotifInOneIntv(veciter(int)& iterStart, veciter(int)& iterEnd,
			i2iHMap& vertex2Pos, DisjointSet*& disjointSet, int& vertexNum,
			vec(int)& combineCCPos, int& realMotifNum,
			i2iHMap& root2Comp, vec(CComponents*)& tempComponents,
			vec(int)& saveCCPos, int motifStartT, int motifEndT,
			vec(TMotif*)*& result, int k, long long& motifNumber, int TFchoice);
		/*generate motifs in one interval for DFTM*/
		void genMotifInOneIntvDynamic(veciter(int)& iterStart, veciter(int)& iterEnd,
			i2iHMap& vertex2Pos, DisjointSet*& disjointSet, int& vertexNum,
			vec(int)& combineCCPos, int& realMotifNum,
			i2iHMap& root2Comp, vec(CComponents*)& tempComponents,
			vec(int)& saveCCPos, int motifStartT, int motifEndT,
			vec(TMotif*)*& result, int k, long long& motifNumber, int TFchoice);

		#pragma region generateMaxTM
		/*fetch new edges from one R edge set and insert into the disjoint set (maintain connectivity)*/
		void maintainUFSet(veciter(int)& infoBegin,
			veciter(int)& infoEnd, i2iHMap& vertex2Pos, DisjointSet*& ufset, int&vertexNum,
			vec(int)& combineCCPos);
		/*combine components which are connected after adding new edges from one R edge set,
		combine components before adding new edges for less computation
		(combine components in generateMaxTM case 3)*/
		void combineComponents(vec(CComponents*)& tempComponents,
			i2iHMap& vertex2Pos, DisjointSet*& disjointSet,
			i2iHMap& root2Comp, /*int& tempComponentsSize,*/int& realMotifNum,
			vec(int)& combineCCPos);
		/*add edge into generated motif or generate new motif and check left expandable
		(generateMaxTM case 1, 2 and add edge into generated motif in generateMaxTM case 3)*/
		void updateNewEdgeInfo(
			veciter(int)& infoBegin, veciter(int)& infoEnd,
			vec(CComponents*)& tempComponents,
			i2iHMap& vertex2Pos, DisjointSet*& disjointSet,
			i2iHMap& root2Comp, /*int& tempComponentsSize,*/ int& realMotifNum,
			vec(int)& saveCCPos, int startTime);
		#pragma endregion

		void generateExpTM(vec(int)&saveCCPos,
			vec(CComponents*)& tempComponents,
			vec(TMotif*)*& result, int k, int motifStartT, int motifEndT,
			long long& motifNumber, int TFchoice);

		void generateExpTMDynamic(vec(int)&saveCCPos,
			vec(CComponents*)& tempComponents,
			vec(TMotif*)*& result, int k, int motifStartT, int motifEndT,
			long long& motifNumber, int TFchoice);
	#pragma endregion 

	#pragma endregion

	/*DFTM (row number<=T-k+1)*/
		void findTMotifsDynamic(int k, vec(TMotif*)*& newResult, int oriEndT,
			i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber, int TFchoice);

	/*computeRES for DFTM (row <= T-k+1)*/
	virtual void computeRESForDFTM(int intvB, int intvE,
		vec(int)*& selectedEdge, int& selectedNum,
		unordered_map<Label, bool>& fixLabel,
		bool isEdgeTypeFixed, int k, vec(TMotif*)*& result, int pos) = 0;


	//reset vertexToId and union-find set
	/*void resetStructure() {
		vertexToId.clear();
		ufset->initialize();
	}*/
protected:
	#pragma region save all information from temporal graph file and create structures(eL table and IC trees)
	/*
	create the eL table and IC trees
	four vectors save four number u,v,t,w respectively
		u,v:two endpoints of an edge
		t: timestamp
		w: edge label
	*/
	virtual void createStructure(vec(int)& u_arr, vec(int)& v_arr,
		vec(int)& t_arr, vec(Label)& w_arr) = 0;

	/*save all information from input temporal graph*/
	void loadInformation(const char*& src, set<Edge>& origEdge,
		vec(int)& u_arr, vec(int)& v_arr,
		vec(int)& t_arr, vec(Label)& w_arr);

	/*save all information from input temporal graph with interval [fixedS, fixedE]*/
	void loadInformation(const char*& src, set<Edge>& origEdge,
		vec(int)& u_arr, vec(int)& v_arr,
		vec(int)& t_arr, vec(Label)& w_arr,
		int fixedS,int fixedE);
	#pragma endregion

	#pragma region auxiliary variable and structure 
		int nTimestamp;	// number of snapshots
		int startT, endT; // beginning timestamp and ending timestamp
		size_t nNode, nEdge; // number of nodes and edges
		Edge* edgeList;// map edges' id to its object
		int maxIntervalLength;//used for reducing memory of R edge set
		DisjointSet* ufset;//disjoint set
		unordered_map<int, int> vertexToId;// map from an node to its id (used for disjoint set)
		map<Edge, int> *edge2ind; // map from an edge to its id (only used for output)

		//used for incremental algorithms to maintain the result 
		i2iHMap* edgeToMotifId; //O(E)
		//int* edgeToMotifId; //O(E)
		//int* minEndTime; //O(max motif number)
		SaveCCInfo* motifSaveInfo; //O(max motif number)

		Intv *EMaxIntvl;
public:
		static int maxMotifNum; 
	#pragma endregion
};
