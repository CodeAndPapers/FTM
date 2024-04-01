#include"DisjointSet.h"
#include"stdafx.h"

#pragma once
/*find the root of the node whose position is num*/
int DisjointSet::find(int num) {
	if (num >= size) return -1;
	if (parent[num] == -1) return -1;
	int root = num;
	while (root != parent[root])
		root = parent[root];
	int temp=num,p;
	while (temp != root) {
		p = parent[temp];
		parent[temp] = root;
		temp = p;
	}
	return parent[num];
}
//
///*union two trees where two nodes locate*/
//void DisjointSet::unionfind(int a, int b) {
//	if (parent[a] == -1 && parent[b] == -1) {
//		parent[a] = parent[b] = b;
//		depth[b] += 1;
//	}
//	else if (parent[a] == -1) {
//		parent[a] = parent[b];
//		depth[b] += 1;
//	}
//	else if (parent[b] == -1) {
//		parent[b] = parent[a];
//		depth[a] += 1;
//	}
//	else {
//		int p_a = find(a);
//		int p_b = find(b);
//		if (p_a == p_b)return;
//		if (depth[p_a] > depth[p_b]) {
//			parent[p_b] = parent[p_a];
//		}
//		else if (depth[p_a] < depth[p_b]) {
//			parent[p_a] = parent[p_b];
//		}
//		else {
//			parent[p_a] = parent[p_b];
//			depth[p_b] += 1;
//		}
//	}
//}

/*union two trees where two nodes locate
	and return the old root of tree which is combined to other tree*/
int DisjointSet::newUnionVertexs(int a, int b) {
	if (parent[a] == -1 && parent[b] == -1) {
		parent[a] = parent[b] = b;
		sizeOfTr[b] += 1;
	}
	else if (parent[a] == -1) {
		int p_b = find(b); 
		parent[a] = parent[p_b];
		sizeOfTr[p_b] += 1;
	}
	else if (parent[b] == -1) {
		int p_a = find(a); 
		parent[b] = parent[p_a];
		sizeOfTr[p_a] += 1;
	}
	else {
		int p_a = find(a);
		int p_b = find(b);
		if (p_a == p_b)return -1;
		int rtn;
		if (sizeOfTr[p_a] > sizeOfTr[p_b]) {
			rtn = parent[p_b];
			parent[p_b] = parent[p_a];
			sizeOfTr[p_a] += sizeOfTr[p_b];
		}
		else if (sizeOfTr[p_a] < sizeOfTr[p_b]) {
			rtn = parent[p_a];
			parent[p_a] = parent[p_b];
			sizeOfTr[p_b] += sizeOfTr[p_a];
		}
		else {
			rtn = parent[p_a];
			parent[p_a] = parent[p_b];
			sizeOfTr[p_b] += sizeOfTr[p_a];
		}
		return rtn;
	}
	return -1;
}