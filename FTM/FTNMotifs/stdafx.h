#pragma once

#include "targetver.h"
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#else
#define DBG_NEW new
#endif

#include <stdio.h>
#include<iostream>
#include<cctype>
#include<cstring>
#include<algorithm>
#include<vector>
#include<queue>
#include<exception>
#include<unordered_map>
#include<string>
#include<fstream>
#include<stack>
#include<map>
#include<set>
#include<ctime>
#include<cstdlib>
#include<Windows.h>
#include<psapi.h>
#include<iomanip>
#include <memory>
#include <cassert>
#include <chrono>
using namespace std;

#pragma region default input parameter
#define DEFAULT_K 10 //default frequent condition k
#define FILE_NAME_LENGTH 100 //maximum length of file name
#define OUTPUT_FILE "example-output.txt" //default output file name
#define INPUT_FILE "example-temporal-graph.txt" //default input temporal graph file name
#define K_FILE "example-k.txt" //default file name for frequent condition k
#define FIXEDLABEL_FILE "example-fixedlabel.txt" //default file name for fixed label
#define ALGORITHM_ID 1 //default algorithm id
#define RUNTIMES 1 //default times of running the algorithm
#define INDEXID 1 //default index id (EL structure)
#pragma endregion

#pragma region util function
#define CLEARALL(a,value,num,type) memset(a,value,sizeof(type)*num)
#define OUTPUT_FILE_OPEN freopen(FindTMotif::outputSrc, "a", stdout); ios::sync_with_stdio(false);
#define OUTPUT_FILE_CLOSE fclose(stdout);
#define VEC_RELEASE(type,vec) vector<type>().swap(vec)
#pragma endregion

#pragma region output parameter
#define EDGE_NUM "edges: "
#define INTV_NUM "maxVe: "
#define MEAN_INTV_NUM "averVe: "
#define MOTIF_ID "No: "
#define MIDRESULT_MEMORY "Test intermediate result memory:\n"
#define RUN_TIMES "run times: "
#define MOTIF_NUM "motifs: "
#define SELECT_EDGE "maxSm: "
#define MEAN_SELECT_EDGE "averSm: "
#pragma endregion


#pragma region error output
#define EXIT exit(-1);
#pragma endregion

#pragma region type
#define vec(a) vector<a>
#define veciter(a) vector<a>::iterator

using Label = int;
struct TEdge {//edge's id and label
	int id;
	Label label;
	TEdge(int i, Label lab) {
		id = i;
		label = lab;
	}
	TEdge() = default;
	bool operator==(const TEdge & tedge) const {
		return (id == tedge.id && label == tedge.label);
	}
	bool operator!=(const TEdge & tedge) const {
		return (id != tedge.id || label != tedge.label);
	}
	bool operator<(const TEdge& tedge) const {
		return id < tedge.id;
	}
	bool operator>(const TEdge& tedge) const {
		return id > tedge.id;
	}
};

struct IntPair {//used as interval[first,second]
	int first;
	int second;
	IntPair(int f, int s) {
		first = f;
		second = s;
	}
	bool operator < (const IntPair & intpair) const {
		return (first < intpair.first || 
			(first == intpair.first)&& second < intpair.second);
	}
	IntPair() {
		first = 0;
		second = 0;
	}
};

using l2lHMap = unordered_map<long long, long long>;
using i2iHMap = unordered_map<int, int>;
using i2bHMap = unordered_map<int, bool>;
using i2iHMap_Iter = unordered_map<int, int>::iterator;
#pragma endregion

//information saved for R edge sets and EMaxIntvl
struct SaveInfoNoStruct {
	int edgeId; //edge in R edge sets
	int startT; //interval for EMaxIntvl[e], only need to save the starting timestamp
	Label label; //label of e in interval
	SaveInfoNoStruct(int id, int start, Label w) {
		edgeId = id;
		startT = start;
		label = w;
	}
};
using SAVEINFONOST = SaveInfoNoStruct;
using SAVEINFONOST_Vec = vec(SAVEINFONOST);
using SAVEINFONOST_VIter = veciter(SAVEINFONOST);

//information saved for R edge sets and EMaxIntvl
struct SaveInfo {
	int edgeId; //edge in R edge sets
	//int startT; //interval for EMaxIntvl[e], only need to save the starting timestamp
	Label label; //label of e in interval
	SaveInfo(int id/*, int start*/, Label w) {
		edgeId = id;
		//startT = start;
		label = w;
	}
};
using SAVEINFO = SaveInfo;
using SAVEINFO_Vec = vec(SAVEINFO);
using SAVEINFO_VIter = veciter(SAVEINFO);