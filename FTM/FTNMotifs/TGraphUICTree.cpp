#include"TGraphUICTree.h"
#include"stdafx.h"

TGraphUICTree::TGraphUICTree(const TGraphUICTree& ances) :TGraph(ances) {
	*tree = *ances.tree;
}

#pragma region construct and update the temporal graph
/* construct the temporal graph
1) graph file: each line of the file describe an edge in a timestamp.
each line has four number: u,v,t,w (separated by ',' with no other space)
for weight w of edge (u,v) in time t. Ids of node u and v are not guaranteed to be
continuous, while the timestamps t are continuous, i.e. all edges in time 0 come
first, and followed by edges in time 1, 2...
*/
void TGraphUICTree::constructGraph(const char* src) {
	vec(int) u_arr, v_arr, t_arr;
	vec(Label) w_arr;
	set<Edge> origEdge;
	loadInformation(src, origEdge, u_arr, v_arr, t_arr, w_arr);
	createStructure(u_arr, v_arr, t_arr, w_arr/*, name2id*/);
}

/*construct the temporal graph, but fixed the startT and endT*/
void TGraphUICTree::constructGraph(const char* src,
	int startT, int endT) {
	vec(int) u_arr, v_arr, t_arr;
	vec(Label) w_arr;
	set<Edge> origEdge;
	loadInformation(src, origEdge, u_arr, v_arr, t_arr, w_arr,
		startT, endT);
	createStructure(u_arr, v_arr, t_arr, w_arr/*, name2id*/);
}

//increase snapshots of the temporal graph
void TGraphUICTree::changeGraph(const char* src, int oldEndT, int limitNewEndT) {
	vec(int) u_arr, v_arr, t_arr;
	vec(Label) w_arr;//allocate the memory
	u_arr.reserve(ALLOC_MEM);
	v_arr.reserve(ALLOC_MEM);
	t_arr.reserve(ALLOC_MEM);
	w_arr.reserve(ALLOC_MEM);
	int u, v, t;
	Label w;

#pragma region load file
	FILE* file;
	file = fopen(src, "r+");
	if (!file) exit(-1);
	char line[LINE_LENGTH];
	CLEARALL(line, 0, LINE_LENGTH, char);
	int sep1, sep2, sep3;//separator pos
	int newEndT = oldEndT, oldTimestamp;
	if (limitNewEndT == -1) limitNewEndT = 0x7fffffff;
	try {
		while (fgets(line, LINE_LENGTH, file)) {
			if (strlen(line) == 0) continue;
			sep1 = (int)(find(line, line + LINE_LENGTH, SEP_CHAR) - line);
			sep2 = (int)(find(line + sep1 + 1, line + LINE_LENGTH, SEP_CHAR) - line);
			sep3 = (int)(find(line + sep2 + 1, line + LINE_LENGTH, SEP_CHAR) - line);
			u = STR2INT(line);
			v = STR2INT(line + sep1 + 1);
			t = STR2INT(line + sep2 + 1);
			w = STR2INT(line + sep3 + 1);

			if (t > oldEndT&&t <= limitNewEndT) {
				if (t > newEndT)newEndT = t;//new endT
				u_arr.emplace_back(u);
				v_arr.emplace_back(v);
				t_arr.emplace_back(t);
				w_arr.emplace_back(w);
			}
		}
		fclose(file);
		this->endT = newEndT;//update endT of graph
		oldTimestamp = this->nTimestamp;
		this->nTimestamp = newEndT - this->startT + 1;
	}
	catch (exception& e) {
		Util::printError(e);
		cout << "load error" << endl;
	}
#pragma endregion
	
	int newTimestamp = newEndT - oldEndT;
	int newStartT = oldEndT + 1;
	Label** eW = DBG_NEW Label*[newTimestamp];//edge's new weight
	for (int i = 0; i < newTimestamp; i++) {
		eW[i] = DBG_NEW Label[nEdge];
		for (int j = 0; j < nEdge; j++) {
			eW[i][j] = 0x7fffffff;
		}
	}
	delete[] EMaxIntvl;
	EMaxIntvl = DBG_NEW Intv[nEdge];

#pragma region prepare for updating
	clock_t startTime, endTime;
	startTime = clock();
	int edgeInd;
	size_t length = u_arr.size();
	for (size_t i = 0; i < length; i++) {
		u = u_arr[i];
		v = v_arr[i];
		t = t_arr[i];
		w = w_arr[i];

		Edge e(u, v);
		edgeInd = edge2ind->find(e)->second;
		if (eW[t - newStartT][edgeInd] == 0x7fffffff)
			eW[t - newStartT][edgeInd] = w;
	}
#pragma endregion

#pragma region update index
	Label nowW, breW;
	int intv_s;//current Intv's start time
	Intv updateIntv;
	Intv intv;
	int intvLen, pos;
	for (int id = 0; id < nEdge; id++) {//O(E)
		tree[id]->search(updateIntv, oldEndT);
		intv_s = updateIntv.startT;
		breW = updateIntv.value;
		for (int i = oldEndT + 1; i <= newEndT; i++) {//O(deltaT)
			pos = i - oldEndT - 1;
			nowW = eW[pos][id];
			if (nowW != breW) {//next Intv
				if (i != oldEndT + 1) {
					intvLen = i - intv_s;
					intv.setValue(intv_s, i - 1, breW);
					maxIntervalLength = max(maxIntervalLength, intvLen);
					if (intv_s <= oldEndT) {//change interval
						tree[id]->updateNode(intv, i - 1);//log(n)
					}
					else {//new interval
						tree[id]->insert(intv);
					}
				}
				intv_s = i;//update
				breW = nowW;
			}
		}
		intvLen = newEndT - intv_s + 1;
		maxIntervalLength = max(maxIntervalLength, intvLen);
		if (intv_s <= oldEndT) {//change interval
			intv.setValue(intv_s, newEndT, breW);
			tree[id]->updateNode(intv, newEndT);//log(n)
		}
		else {//new interval
			intv.setValue(intv_s, newEndT, breW);
			tree[id]->insert(intv);
		}
	}
	endTime = clock();
	cout << "update: " << endTime - startTime << "ms" << endl;

	for (int i = 0; i < newTimestamp; i++) {
		delete[]eW[i];
	}
	delete[]eW;
#pragma endregion
}
#pragma endregion

/*
create the eL table and IC trees
four vectors save four number u,v,t,w respectively
	u,v:two endpoints of an edge
	t: timestamp
	w: edge label
*/
void TGraphUICTree::createStructure(vec(int)& u_arr, vec(int)& v_arr,
	vec(int)& t_arr, vec(Label)& w_arr) {
	Label** eW = DBG_NEW Label*[nTimestamp];
	for (int i = 0; i < nTimestamp; i++) {
		eW[i] = DBG_NEW Label[nEdge];
		for (int j = 0; j < nEdge; j++) {
			eW[i][j] = 0x7fffffff;
		}
	}
	EMaxIntvl = DBG_NEW Intv[nEdge];
	for (int j = 0; j < nEdge; j++) {
		EMaxIntvl[j].endT = EMaxIntvl[j].startT = -1;
	}
	tree = DBG_NEW ICTree*[nEdge];//use Intv tree

	clock_t startTime, endTime;
	startTime = clock();
	int u, v, t, edgeInd;
	Label w;
	size_t length = u_arr.size();
	for (size_t i = 0; i < length; i++) {
		u = u_arr[i];
		v = v_arr[i];
		t = t_arr[i];
		w = w_arr[i];

		Edge e1(u, v);
		edgeInd = edge2ind->find(e1)->second;
		if (eW[t - startT][edgeInd] == 0x7fffffff)
			eW[t - startT][edgeInd] = w;
	}

	Label nowW, breW;
	int intervalLen;
	int intv_s;//current Intv's start time
	//oriSearchT = 0;
	int intvNum = 0;// number of intervals
	long long intvSum = 0;// sum of intervals number
	int maxIntvNum = 0;// maximum number of intervals
	for (int id = 0; id < nEdge; id++) {//O(E)
		breW = eW[0][id];
		intv_s = startT;
		tree[id] = DBG_NEW ICTree();//create the tree 
		for (int i = startT + 1; i <= endT; i++) {//O(T)
			nowW = eW[i - startT][id];
			if (nowW != breW) {//next Intv
				tree[id]->insert(Intv(intv_s, i - 1, breW));
				intvNum++;
				intervalLen = i - intv_s;
				maxIntervalLength = max(maxIntervalLength, intervalLen);
				intv_s = i;//update
				breW = nowW;
			}
		}
		tree[id]->insert(Intv(intv_s, endT, breW));
		intvNum++;
		intvSum += intvNum;
		maxIntvNum = max(maxIntvNum, intvNum);
		intvNum = 0;
		intervalLen = endT - intv_s + 1;
		maxIntervalLength = max(maxIntervalLength, intervalLen);
	}
	cout << INTV_NUM << maxIntvNum << endl;
	cout << MEAN_INTV_NUM << intvSum * 1.0 / nEdge << endl;
	cout << "Balance trees: " << (37 * intvSum + 8 * (nEdge + 1)) << "Byte" << endl;
	endTime = clock();
	cout << "load: " << endTime - startTime << "ms" << endl;

	for (int i = 0; i < nTimestamp; i++) {
		delete[]eW[i];
	}
	delete[]eW;
}

#pragma region key static algorithms
/*FTM*/
void TGraphUICTree::findTMotifs(int k, vec(TMotif*)*& result,
	i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber,
	int choiceStartT, int choiceEndT) {
#pragma region initialize
	int lastTime = choiceEndT - k + 1;//last start time 
	i2iHMap vertex2Pos;//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp;

	vec(CComponents*) tempComponents;//temporary component list CC
	vec(int)* edgeSetsR = DBG_NEW vec(int)[maxIntervalLength];// R edge sets
	DisjointSet* disjointSet;//disjoint set

	veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
	veciter(CComponents*) listEnd;//iterator for tempComponents
	int selectedNum;//the number of edges in edgeSetsR

	vec(int) saveMotifPos;//record the position of components which need to be saved
	/*record the position of components which need to combine with other components*/
	vec(int) combineMotifPos;

	long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
	long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
	clock_t begin;//record time
	long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]

	int vertexNum;//the number of vertexs
	int realMotifNum;//the number of generated connected component
	int edgeEndT;//record the currently considering end time of R edge sets
	int Te;//the end time of interval
#pragma endregion

	for (int Ts = choiceStartT; Ts <= lastTime; Ts++) {//O(T)
		// select edges which keep fixed in [Ts,Ts+k-1]
		begin = clock();
		selectedNum = 0;
		Te = Ts + k - 1;
		computeRES(Ts, Te,
			edgeSetsR, selectedNum,
			fixLabel, isEdgeTypeFixed);

		//initalize for every row
		maxSelect = max(maxSelect, selectedNum);
		selectedSum += selectedNum;
		disjointSet = DBG_NEW DisjointSet((int)((selectedNum << 1) + 1));
		vertex2Pos.clear();
		root2Comp.clear();
		vertexNum = 0;
		realMotifNum = 0;
		edgeEndT = choiceEndT;
		Test::compr += clock() - begin;

		//Test S(m,T)
		if (maxIntervalLength <= choiceEndT - Te) smt = 0;
		else smt = edgeSetsR[choiceEndT - Te].size();
		maxSmt = max(maxSmt, smt);
		allSmt += smt;

		//R2L
		for (int tempT = choiceEndT - Te; tempT >= 0; tempT--, edgeEndT--) {
			if (tempT >= maxIntervalLength) continue;
			iterStart = edgeSetsR[tempT].begin();
			iterEnd = edgeSetsR[tempT].end();
			if (iterStart == iterEnd) continue;
			//generateMaxTM and generatedExpTM
			genMotifInOneIntv(iterStart, iterEnd,
				vertex2Pos, disjointSet, vertexNum,
				combineMotifPos, realMotifNum, root2Comp, tempComponents,
				saveMotifPos, Ts, edgeEndT,
				result, k, motifNumber);
			edgeSetsR[tempT].clear();
		}

		//testing
		Test::updateMemoryUse();
		//Test::showMemoryUse();

		//release
		delete disjointSet;
		listEnd = tempComponents.end();
		for (veciter(CComponents*) listIter = tempComponents.begin();
			listIter != listEnd; ++listIter) {
			if (*listIter != NULL)
				delete *listIter;
		}
		tempComponents.clear();
	}
	delete[] edgeSetsR;
	cout << SELECT_EDGE << maxSelect << endl;
	cout << MEAN_SELECT_EDGE << selectedSum / (lastTime - choiceStartT + 1) << endl;
	cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime - choiceStartT + 1) << endl;
}
#pragma endregion


/*compRES for FTM*/
void TGraphUICTree::computeRES(int intvB, int intvE,
	vec(int)*& edgeSetsR, int& selectedNum,
	unordered_map<int, bool>& fixLabel, bool isEdgeTypeFixed) {
	//int timePos = intvB - startT;
	Label edgeType;
	//unordered_map<int, bool>::iterator fixLabelEnd = fixLabel.end();
	int intvLen = intvE - intvB + 1;
	int rtnStartT, rtnEndT;

	for (int i = 0; i < nEdge; i++) {
		if (EMaxIntvl[i].endT >= intvE && EMaxIntvl[i].startT <= intvB) {//case 4
			edgeType = EMaxIntvl[i].value;
			/*if (!isEdgeTypeFixed ||
				fixLabel.find(edgeType) != fixLabelEnd) {*/
				selectedNum++;
				edgeSetsR[EMaxIntvl[i].endT - intvE].
					emplace_back(i);
			/*}*/
		}
		else if (EMaxIntvl[i].endT >= intvB && EMaxIntvl[i].startT <= intvE) {//case 3
			continue;
		}
		else {//case 1 and 2
			bool flag = tree[i]->containQuery(intvB, intvE, edgeType, rtnStartT, rtnEndT);
			EMaxIntvl[i].endT = rtnEndT;
			EMaxIntvl[i].startT = rtnStartT;
			EMaxIntvl[i].value = edgeType;
			if (flag && rtnEndT >= intvE/*&& (!isEdgeTypeFixed ||
					fixLabel.find(edgeType) != fixLabelEnd)*/) {
				selectedNum++;
				edgeSetsR[rtnEndT - intvE].emplace_back(i);
			}
		}
	}
}


/*get the label of edge with edgeId at the time*/
Label TGraphUICTree::getWeight(int time, int edgeId) {
	Intv intv;
	tree[edgeId]->search(intv, time);
	return intv.value;
}

/*compRES for DFTM (row number<=T-k+1)*/
void TGraphUICTree::computeRESForDFTM(int intvB, int intvE,
	vec(int)*& edgeSetsR, int& selectedNum,
	unordered_map<Label, bool>& fixLabel,
	bool isEdgeTypeFixed, int k, vec(TMotif*)*& allResult, int pos) {
	vec(TMotif*)& result = allResult[pos];
	//int timePos = intvE - startT;
	Label edgeType;
	//unordered_map<int, bool>::iterator fixLabelEnd = fixLabel.end();
	int intvLen = intvE - intvB + 1;
	int endTime = intvB + k - 1;
	int rtnStartT, rtnEndT;
	int check;
	veciter(TMotif*) iterEnd = result.end();
	TMotif* tempMotif;
	vec(TEdge)* motifEdge; 
	SaveCCInfo saveinfo;
	vec(SaveCCInfo)* saveinfoList;
	veciter(TEdge) edgeIterEnd;
	int id;
	int minEndT, minEndTEdge;
	int motifPos = 0;
	for (auto iter = result.begin(); iter != iterEnd; ++iter, ++motifPos) {//consider edges only in TF^o[m,T]
		tempMotif = *iter;
		motifEdge = tempMotif->getMotifEdge();
		edgeIterEnd = motifEdge->end();
		minEndT = 0x7fffffff;
		
		saveinfoList = tempMotif->getOtherEdge();
		if (saveinfoList != NULL) saveinfo = *saveinfoList->begin();
		else { saveinfo.saveEndTime = -1; }
		motifSaveInfo[motifPos] = saveinfo;
		for (auto edgeIter = motifEdge->begin();
			edgeIter != edgeIterEnd; ++edgeIter) {
			id = edgeIter->id;
			//edgeToMotifId[id] = motifPos;

			if (EMaxIntvl[id].endT >= intvE && EMaxIntvl[id].startT <= intvB) {//case 4
				edgeType = EMaxIntvl[id].value;
				if (/*(!isEdgeTypeFixed ||
					fixLabel.find(edgeType) != fixLabelEnd)
					&&*/ EMaxIntvl[id].endT >= intvE) {
					check = EMaxIntvl[id].endT - intvE;
					selectedNum++;
					if (minEndT >= check) {
						minEndT = check;
						minEndTEdge = id;
					}
					edgeSetsR[check].
						emplace_back(id);
				}
			}
			else if (EMaxIntvl[id].endT >= intvB && EMaxIntvl[id].startT <= intvE) {//case 3
				continue;
			}
			else {//case 1 and 2
				bool flag = tree[id]->containQuery(intvB, endTime, edgeType, rtnStartT, rtnEndT);
				EMaxIntvl[id].endT = rtnEndT;
				EMaxIntvl[id].startT = rtnStartT;
				EMaxIntvl[id].value = edgeType;
				if (flag/* && (!isEdgeTypeFixed ||
						fixLabel.find(edgeType) != fixLabelEnd)*/) {
					if (rtnEndT >= intvE) {
						selectedNum++; 
						check = rtnEndT - intvE;
						if (minEndT >= check) {
							minEndT = check;
							minEndTEdge = id;
						}
						edgeSetsR[check].emplace_back(id);
					}
				}
			}
		}
		//for (auto edgeIter = motifEdge->begin();
		//	edgeIter != edgeIterEnd; ++edgeIter) {
		//	id = edgeIter->id;
		//minEndTime[motifPos] = minEndT + intvE;
		//}
		edgeToMotifId->insert_or_assign(minEndTEdge, motifPos);
		delete tempMotif;
	}
}
