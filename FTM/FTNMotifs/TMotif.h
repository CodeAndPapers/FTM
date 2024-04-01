#pragma once
#include "stdafx.h"
#include "DisjointSet.h"
#include "Edge.h"


struct SaveCCInfo {
	int savePos;//cc's saving position in the result
	int saveEndTime;//the end time of the interval of cc
	int saveCCInfoPos;//cc's saveCCInfo position
	SaveCCInfo(int saveP, int saveET) {
		savePos = saveP;
		saveEndTime = saveET;
		saveCCInfoPos = -1;
	}
	SaveCCInfo(int saveP, int saveET, int saveCCInfoP) {
		savePos = saveP;
		saveEndTime = saveET;
		saveCCInfoPos = saveCCInfoP;
	}
	
	SaveCCInfo() {
		savePos = saveEndTime = saveCCInfoPos = -1;
	}
};

/*temporal motifs*/
class TMotif {
public:
	TMotif():startT(0),endT(0),motifEdge(NULL), otherEdge(NULL) {}

	TMotif(vec(TEdge)& selectedEdge,int start,int end) {
		startT = start;
		endT = end;
		motifEdge = NULL;
		otherEdge = NULL;
		setMotifEdge(selectedEdge);
	}
	
	TMotif(int start,int end):startT(start),endT(end), otherEdge(NULL)
	{
		motifEdge = DBG_NEW vec(TEdge)();
	}

	// copy-constructed
	TMotif(const TMotif& motif) {
		this->motifEdge = DBG_NEW vec(TEdge)(*motif.motifEdge);
		this->otherEdge = DBG_NEW vec(SaveCCInfo)(*motif.otherEdge);
		startT = motif.startT;
		endT = motif.endT;
		edgeNumber = motif.edgeNumber;
	}

	/* add edge to the motif*/
	inline void addEdge(int id, Label edgeW) {
		motifEdge->emplace_back(id, edgeW);
	}

	inline void updateEdgeNumber() {
		edgeNumber = motifEdge->size();
	}

	//size
	//inline size_t getSize() {
	//	return motifEdge->size();
	//}

	// /* transform from CComponents to TMotif 
	//	and update the ufset*/
	//void transform(CComponents*& tempMotif, Edge*& edgeList,
	//	i2iHMap& vertex2Pos, DisjointSet*& ufset, int& vertexNum) {
	//	veciter(TEdge) edgeEnd = motifEdge->end();
	//	
	//	tempMotif = DBG_NEW CComponents();
	//	Edge* temp;
	//	int motifPos;
	//	int sId, tId, vertex;//vertexs'id in the union-find set
	//	i2iHMap_Iter vertexIt;
	//	for (auto iter = motifEdge->begin();
	//		iter != edgeEnd; ++iter) {
	//		tempMotif->edges.emplace_back(*iter);
	//		temp = &edgeList[iter->id];//new inserted edge
	//		vertex = temp->s;
	//		vertexIt = vertex2Pos.find(vertex);//O(1)
	//		if (vertexIt == vertex2Pos.end()) {
	//			sId = vertexNum;
	//			vertex2Pos[vertex] = vertexNum++;
	//		}
	//		else sId = vertexIt->second;
	//		vertex = temp->t;
	//		vertexIt = vertex2Pos.find(vertex);//O(1)
	//		if (vertexIt == vertex2Pos.end()) {
	//			tId = vertexNum;
	//			vertex2Pos[vertex] = vertexNum++;
	//		}
	//		else tId = vertexIt->second;
	//		
	//		/*union-find operation means that sId and tId are connected*/
	//		motifPos = ufset->newUnionVertexs(sId, tId);
	//	}
	//	tempMotif->root = ufset->find(sId);
	//	tempMotif->startT = this->startT;
	//}

	// /* output the information of TMotif*/
	//friend ostream& operator<<(ostream& output, const TMotif& motif) {
	//	//output << "output temporal motif information:" << endl;
	//	vector<TEdge>* vec = motif.motifEdge;
	//	output << "startT: " << motif.startT << "\tendT: " << motif.endT << endl;
	//	output << EDGE_NUM << vec->size();
	//	
	//	vector<TEdge>::iterator it = vec->begin();
	//	for (; it != vec->end(); ++it) {
	//		output << "\n"<<*it;
	//	}
	//	return output;
	//}

	vec(TEdge)* getMotifEdge() { return motifEdge; }

	vec(SaveCCInfo)* getOtherEdge() { return otherEdge; }

	inline int getStartT() const { return startT; }

	inline int getEndT() const { return endT; }

	//size
	inline size_t getSize() {
		return edgeNumber;
	}

	void setMotifEdge(vec(TEdge)& motifEdge) {
		if(this->motifEdge)
			delete this->motifEdge;
		this->motifEdge = DBG_NEW vec(TEdge)(motifEdge);
		this->edgeNumber = this->motifEdge->size();
	}

	~TMotif() {
		if(motifEdge){
			motifEdge->clear();
			delete motifEdge;
		}
		if (otherEdge) {
			otherEdge->clear();
			delete otherEdge;
		}
	}
	
	/* linked to other motifs*/
	inline int linkToMotifs(SaveCCInfo& saveInfo, TMotif*& motif) {
		if (this->otherEdge == NULL)
			this->otherEdge = DBG_NEW vec(SaveCCInfo)();
		otherEdge->emplace_back(saveInfo);
		this->edgeNumber += motif->getSize();
		return otherEdge->size() - 1;
	}


	/* linked to other motifs for motifs with interval = [.,endT], in order to update SaveCCInfo of those motifs when they (interval = [.,endT]) update in the incremental algorithm*/
	inline void tempLinkToMotifs(SaveCCInfo& saveInfo) {
		if (this->otherEdge == NULL)
			this->otherEdge = DBG_NEW vec(SaveCCInfo)();
		otherEdge->emplace_back(saveInfo);
	}

	inline void setInterval(int startT, int endT) {
		this->startT = startT;
		this->endT = endT;
	}

	inline void setEndT(int endT) {
		this->endT = endT;
	}

	friend bool operator==(const TMotif&a,const TMotif&b) {
		if (a.startT == b.startT&&a.endT == b.endT) {
			vec(TEdge)* aEdge = a.motifEdge,* bEdge=b.motifEdge;
			bool same = true;
			for (auto iter = aEdge->begin();
				iter != aEdge->end(); ++iter) {
				if (find(bEdge->begin(), bEdge->end(), *iter) == bEdge->end()) {
					same = false;
					break;
				}
			}
			return same;
		}
		return false;
	}

private:
	vec(TEdge)* motifEdge; //motif's edges with labels
	int startT, endT; // starting timestamp(include) and ending timestamp(include)
	int edgeNumber;
	//for memory compression
	vec(SaveCCInfo)* otherEdge; //reuse other motif's edges (already inserted)
};

/*connected components*/
class CComponents {
public:
	vec(TEdge) edges;//cc's edges with labels
	vec(SaveCCInfo) saveInfo;//cc's saving information in the result
	int startT;//the starting timestamp of cc.intvl
	int root;//the root of cc in disjoint set
	CComponents(int start, int root) :startT(start),
		root(root) {}
	CComponents() = default;
	~CComponents() = default;
};