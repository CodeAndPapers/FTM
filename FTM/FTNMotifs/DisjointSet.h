#pragma once
#include "stdafx.h"
/*disjoint set used for maintaining connected components
	(use path compression and union by rank)*/

class DisjointSet {
private:
		int* parent;
		int* sizeOfTr;//the depth of subtree
		int size;
public:
	DisjointSet(int size)
		:size(size){
		parent = DBG_NEW int[size];
		sizeOfTr = DBG_NEW int[size];
		initialize();
	}
	~DisjointSet() {
		delete[]parent;
		delete[]sizeOfTr;
	}

	/*find the root of the node whose position is num*/
	int find(int num);
	
	/*union two trees where two nodes locate*/
	//void unionfind(int a, int b);
	
	/*union two trees where two nodes locate
	and return the old root of tree which is combined to other tree*/
	int newUnionVertexs(int a, int b);

	inline int getSize()const { return size; }

	int* getParent()const { return parent; }

	int* getDepth()const { return sizeOfTr; }

	void initialize() {
		for (int i = 0; i < size; i++) {
			parent[i] = -1;
			sizeOfTr[i] = 1;
		}
	}

	void clear() {
		for (int i = 0; i < size; i++) {
			parent[i] = -1;
			sizeOfTr[i] = 1;
		}
	}

	DisjointSet(const DisjointSet& ufset) {
		size = ufset.getSize();
		parent = DBG_NEW int[size];
		sizeOfTr = DBG_NEW int[size];
		int* d = ufset.getDepth();
		int* p = ufset.getParent();
		for (int i = 0; i < size; ++i) {
			sizeOfTr[i] = d[i];
			parent[i] = p[i];
		}
	}

	void print() {
		for (int i = 0; i < size; i++) {
			cout<<i<<"(p:"<<parent[i]<<
				",d:"<<sizeOfTr[i]<<")"<<endl;
		}
	}
};