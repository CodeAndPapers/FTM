#include "TGraphUEL.h"
#include "stdafx.h"

TGraphUEL::TGraphUEL(const TGraphUEL& ances):TGraph(ances) {
	eW = DBG_NEW Label*[nTimestamp];
	eL = DBG_NEW int*[nTimestamp];
	//eminL = new int[nTimestamp];
	for (int i = 0; i < nTimestamp; i++) {
		eW[i] = DBG_NEW Label[nEdge];
		eL[i] = DBG_NEW int[nEdge];
		//eminL[i] = ances.eminL[i];
		for (int j = 0; j < nEdge; j++) {
			eW[i][j] = ances.eW[i][j];
			eL[i][j] = ances.eL[i][j];
		}
	}
}

#pragma region construct and update the temporal graph
/* construct the temporal graph
1) graph file src: each line of the file describe an edge in a timestamp.
each line has four number: u,v,t,w (separated by ',' with no other space)
for weight w of edge (u,v) in time t. Ids of node u and v are not guaranteed to be
continuous, while the timestamps t are continuous, i.e. all edges in time 0 come
first, and followed by edges in time 1, 2...
*/
void TGraphUEL::constructGraph(const char* src) {
	vec(int) u_arr, v_arr, t_arr;
	vec(Label) w_arr;
	set<Edge> origEdge;
	loadInformation(src, origEdge, u_arr, v_arr, t_arr, w_arr);
	createStructure(u_arr, v_arr, t_arr, w_arr/*, name2id*/);
}

/*construct the temporal graph, but fixed the startT and endT*/
void TGraphUEL::constructGraph(const char* src,
	int startT, int endT) {
	vec(int) u_arr, v_arr, t_arr;
	vec(Label) w_arr;
	set<Edge> origEdge;
	loadInformation(src, origEdge, u_arr, v_arr, t_arr, w_arr,
		startT, endT);
	createStructure(u_arr, v_arr, t_arr, w_arr/*, name2id*/);
	
}

//increase snapshots of the temporal graph
void TGraphUEL::changeGraph(const char* src, int oldEndT, int limitNewEndT) {
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

			if (t > oldEndT && t <= limitNewEndT) {
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

#pragma region prepare for updating
	
	Label** neweW = DBG_NEW Label*[nTimestamp];//edge's new weight
	int** neweL = DBG_NEW int*[nTimestamp];
	int* newMaxeL = DBG_NEW int[nTimestamp];
	CLEARALL(newMaxeL, 0, this->nTimestamp, int);
	for (int i = 0; i < oldTimestamp; i++) newMaxeL[i] = maxeL[i];
	delete[] maxeL;
	maxeL = newMaxeL;
	delete[] EMaxIntvl;
	EMaxIntvl = DBG_NEW Intv[nEdge];
	for (int j = 0; j < nEdge; j++) {
		EMaxIntvl[j].endT = EMaxIntvl[j].startT = -1;
	}
	for (int i = 0; i < nTimestamp; i++) {
		neweW[i] = DBG_NEW Label[nEdge]; 
		neweL[i] = DBG_NEW int[nEdge];
		if (i < oldTimestamp) {
			for (int j = 0; j < nEdge; j++) {
				neweW[i][j] = eW[i][j]; 
				neweL[i][j] = eL[i][j];
			}
			delete[] eW[i];
			delete[] eL[i];
		}
		else {
			CLEARALL(neweL[i], 0, nEdge, int);
			for (int j = 0; j < nEdge; j++) {
				neweW[i][j] = 0x7fffffff;
			}
		}
	}
	delete[] eW;
	delete[] eL;
	eW = neweW;
	eL = neweL;
#pragma endregion

#pragma region update index
	clock_t startTime, endTime;
	startTime = clock();
	int edgeInd, pos;
	size_t length = u_arr.size();
	for (size_t i = 0; i < length; i++) {
		u = u_arr[i];
		v = v_arr[i];
		t = t_arr[i];
		w = w_arr[i];

		Edge e(u, v);
		pos = t - startT;
		edgeInd = edge2ind->find(e)->second;
		if (eW[pos][edgeInd] == 0x7fffffff) {
			eW[pos][edgeInd] = w;
			if (w == eW[pos - 1][edgeInd]) {
				if (eL[pos - 1][edgeInd] == -1) {
					eL[pos - 1][edgeInd] = 1;
				}
				eL[pos][edgeInd] = eL[pos - 1][edgeInd] + 1;
				maxeL[pos] = max(eL[pos][edgeInd], maxeL[pos]);
				maxIntervalLength = max(maxIntervalLength, maxeL[pos]);

				if (t == endT) {
					int first = pos - eL[pos][edgeInd] + 1;
					eL[first][edgeInd] = -eL[pos][edgeInd];
				}
			}
			else {
				int first = pos - 1 - eL[pos - 1][edgeInd] + 1;
				eL[first][edgeInd] = -eL[pos - 1][edgeInd];

				if (t == endT) {
					eL[pos][edgeInd] = -1;
				}
				else eL[pos][edgeInd] = 1;
			}
		}
	}
	
	endTime = clock(); 
	cout << "update: " << endTime - startTime << "ms" << endl;
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
void TGraphUEL::createStructure(vec(int)& u_arr, vec(int)& v_arr,
	vec(int)& t_arr, vec(Label)& w_arr) {
	eW = DBG_NEW Label*[nTimestamp];
	eL = DBG_NEW int*[nTimestamp];
	maxeL = DBG_NEW int[nTimestamp];
	EMaxIntvl = DBG_NEW Intv[nEdge];
	for (int i = 0; i < nTimestamp; i++) {
		eW[i] = DBG_NEW Label[nEdge];
		eL[i] = DBG_NEW int[nEdge];
		maxeL[i] = 1;
		for (int j = 0; j < nEdge; j++) {
			eW[i][j] = 0x7fffffff;
		}
	}
	for (int j = 0; j < nEdge; j++) {
		EMaxIntvl[j].endT = EMaxIntvl[j].startT = -1;
	}

	clock_t startTime, endTime;
	startTime = clock();
	int u, v, t, edgeInd, pos;
	Label w;
	size_t length = u_arr.size();
	//int id_u, id_v;
	for (size_t i = 0; i < length; i++) {//O(TE)
		u = u_arr[i];
		v = v_arr[i];
		t = t_arr[i];
		w = w_arr[i];
		//id_u = name2id->find(u)->second;
		//id_v = name2id->find(v)->second;
		//id2name[id_u] = u;
		//id2name[id_v] = v;
		Edge e1(u, v);
		edgeInd = edge2ind->find(e1)->second;

		pos = t - startT;
		if (eW[pos][edgeInd] == 0x7fffffff) {
			eW[pos][edgeInd] = w;
			if (t == startT) {
				eL[0][edgeInd] = 1;
			}
			else if (w == eW[pos - 1][edgeInd]) {
				eL[pos][edgeInd] = eL[pos - 1][edgeInd] + 1;
				maxeL[pos] = max(eL[pos][edgeInd], maxeL[pos]);
				maxIntervalLength = max(maxIntervalLength, maxeL[pos]);
				if (t == endT) {
					int first = pos - eL[pos][edgeInd] + 1;
					eL[first][edgeInd] = -eL[pos][edgeInd];
				}
			}
			else {
				int first = pos - 1 - eL[pos-1][edgeInd] + 1;
				eL[first][edgeInd] = -eL[pos - 1][edgeInd];
				if (t == endT) {
					eL[pos][edgeInd] = -1;
				}
				else {
					eL[pos][edgeInd] = 1;
				}
			}
		}
	}
	cout << "EL structures: " << 8 * nTimestamp * nEdge << "Byte" << endl;
	endTime = clock();
	cout << "load: " << endTime - startTime << "ms" << endl;
}


#pragma region key static algorithms
/*FTM*/
void TGraphUEL::findTMotifs(int k, vec(TMotif*)*& result,
	i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber,
	int choiceStartT, int choiceEndT, int TFchoice) {
#pragma region initialize
	int lastTime = choiceEndT - k + 1;//last start time 
	i2iHMap vertex2Pos;//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp;

	vec(CComponents*) tempComponents;//temporary component list CC
	int lineNum = lastTime - choiceStartT + 1;
	//int setsRNum = (lineNum * (lineNum + 1)) >> 1;
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
	//int tempPosForEdgeSetsR = 0;
#pragma endregion
	for (int Ts = choiceStartT; Ts <= lastTime; Ts++/*, tempPosForEdgeSetsR += lineNum, lineNum--*/) {//O(T)

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
		//smt = edgeSetsR[choiceEndT - Te].size();
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
				result, k, motifNumber, TFchoice);
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

/*add edge into generated motif or generate new motif and check left expandable
		(generateMaxTM case 1, 2 and add edge into generated motif in generateMaxTM case 3)*/
void TGraphUEL::updateNewEdgeInfo(
	veciter(int)& infoBegin, veciter(int)& infoEnd,
	vec(CComponents*)& tempComponents,
	i2iHMap& vertex2Pos, DisjointSet*& disjointSet,
	i2iHMap& root2Comp, /*int& tempComponentsSize,*/ int& realMotifNum,
	vec(int)& saveCCPos, int startTime) {

#pragma region initialize
	i2bHMap hasSaved;//means whether the motif has been saved 
	i2iHMap_Iter mapIter;//root2Comp's iterator 
	int root;//the root of one vertex in the disjoint set
	int id;//the edge e's id
	int checkStartT;//the edge e's start time of EMaxIntvl[e]
	int label;//the edge e's label
	int ccPos;//cc's position in tempComponents
	/*the cc which has already generated/will generate*/
	CComponents* generatedCC, *newCC;
	int tempComponentsSize;//size of tempComponents
#pragma endregion

	for (veciter(int) infoIter = infoBegin;
		infoIter != infoEnd; ++infoIter) {//new edges in R edge sets

		id = *infoIter;//edge's id

		/*the root of the edge's vertex in the disjoint set*/
		root = disjointSet->find(vertex2Pos[edgeList[id].s]);

		mapIter = root2Comp.find(root);//check which cc the edge belongs to 
		checkStartT = EMaxIntvl[id].startT;//start time of the edge
		//startT = infoIter->startT;//start time of the edge
		label = EMaxIntvl[id].value; //label of the edge

		if (mapIter == root2Comp.end()) {//generateMaxTM case 1
			newCC = DBG_NEW CComponents(checkStartT, root);
			//newCC->edges = DBG_NEW vector<TEdge>();
			newCC->edges.emplace_back(id, label);
			tempComponentsSize = (int)tempComponents.size();
			tempComponents.emplace_back(newCC);
			root2Comp[root] = tempComponentsSize;//update root2Comp

			if (checkStartT == startTime) {//check left unexpandable
				saveCCPos.emplace_back(tempComponentsSize);
				hasSaved[tempComponentsSize] = true;
			}

			realMotifNum++;
		}
		else {//generateMaxTM case 2
			ccPos = mapIter->second;
			generatedCC = tempComponents[ccPos];
			generatedCC->edges.emplace_back(id, label);
			if (generatedCC->startT < checkStartT) {
				generatedCC->startT = checkStartT;
			}
			/* check whether this cc has already been inserted into save list */
			if (hasSaved.find(ccPos) == hasSaved.end()) {
				if (generatedCC->startT == startTime) {//check left unexpandable
					saveCCPos.emplace_back(ccPos);
					hasSaved[ccPos] = true;
				}
			}
		}
	}
}

/*compRES for FTM*/
void TGraphUEL::computeRES(int intvB, int intvE,
	vec(int)*& edgeSetsR, int& selectedNum,
	unordered_map<int, bool>& fixLabel, bool isEdgeTypeFixed) {
	int timePos = intvE - startT;
	int beginPos = intvB - startT, realBeginPos;
	int check;
	Label edgeType;
	unordered_map<int, bool>::iterator fixLabelEnd = fixLabel.end();
	int intvLen = intvE - intvB + 1;
	int lastPos = endT - startT;
	//int intvStart;
	//if (maxeL[timePos] < intvLen) return;
	for (int i = 0; i < nEdge; i++) {
		/*if (isEdgeTypeFixed){
			edgeType = eW[beginPos][i];
			if(fixLabel.find(edgeType) == fixLabelEnd) 
				continue;
		}*/
		
		if (EMaxIntvl[i].endT >= timePos && EMaxIntvl[i].startT <= beginPos) {//case 4
			check = EMaxIntvl[i].endT - timePos;
			edgeSetsR[check].emplace_back(i);
			selectedNum++;
		}
		else if (EMaxIntvl[i].endT >= beginPos && EMaxIntvl[i].startT <= timePos) {//case 3
			continue;
		}
		else {//case 1 and 2
			/*not exists the maximum interval containing [intvB,intvE] for case 1 and 2
			put here for less time*/
			int eLvalue = -eL[beginPos][i];
			realBeginPos = beginPos;
			if (eLvalue < 0) {//not the start time of the interval in incremental algorithm
				realBeginPos = beginPos + eLvalue + 1;
				eLvalue = -eL[realBeginPos][i];
			}
			
			EMaxIntvl[i].startT = realBeginPos;
			EMaxIntvl[i].endT = realBeginPos + eLvalue - 1;

			int check = EMaxIntvl[i].endT - timePos;
			if (check < 0) {
				continue;
			}
			EMaxIntvl[i].value = eW[realBeginPos][i];
			selectedNum++;
			edgeSetsR[check].emplace_back(i);
		}
	}
}



/*get the label of edge with edgeId at the time*/
Label TGraphUEL::getWeight(int time, int edgeId) {
	return eW[time - startT][edgeId];
}


/*compRES for DFTM (row number<=T-k+1)*/
void TGraphUEL::computeRESForDFTM(int intvB, int intvE,
	vec(int)*& edgeSetsR, int& selectedNum,
	 unordered_map<Label, bool>& fixLabel,
	bool isEdgeTypeFixed, int k, vec(TMotif*)*& allResult,int pos) {
	vec(TMotif*)& result = allResult[pos];
	int timePos = intvE - startT;
	int beginPos = intvB - startT;
	//int timeLen = endT - intvE + 1;
	int check;
	unordered_map<int, bool>::iterator fixLabelEnd = fixLabel.end();
	int intvLen = intvE - intvB + 1;
	/*int lastPos = endT - startT;
	int maxIter=0,intvStart;*/
	//if (maxeL[timePos] < intvLen) return;
	veciter(TMotif*) iterEnd = result.end();
	TMotif* tempMotif;
	vec(TEdge)* motifEdge;
	SaveCCInfo saveinfo;
	vec(SaveCCInfo)* saveinfoList;
	veciter(TEdge) edgeIterEnd;
	int id;
	int minEndT, minEndTEdge;
	int motifPos = 0;
	for (auto iter = result.begin();
		iter != iterEnd; ++iter,++motifPos) {//consider edges only in TF^o[m,T]
		tempMotif = *iter;
		motifEdge = tempMotif->getMotifEdge();
		edgeIterEnd = motifEdge->end();
		minEndT = 0x7fffffff;
		
		saveinfoList = tempMotif->getOtherEdge();
		if (saveinfoList != NULL) saveinfo = *saveinfoList->begin();
		else { saveinfo.saveEndTime = -1;}
		motifSaveInfo[motifPos] = saveinfo;
		
		for (auto edgeIter = motifEdge->begin();
			edgeIter != edgeIterEnd; ++edgeIter) {
			id = edgeIter->id;
			//edgeToMotifId[id] = motifPos;
			
			/*not exists the maximum interval containing [intvB,intvE] for case 1 and 2
				put here for less time*/
			if (EMaxIntvl[id].endT >= timePos && EMaxIntvl[id].startT <= beginPos) {//case 4
				check = EMaxIntvl[id].endT - timePos;
				selectedNum++;
				if (minEndT >= check) {
					minEndT = check;
					minEndTEdge = id;
				}
				edgeSetsR[check].emplace_back(id);
			}
			else if (EMaxIntvl[id].endT >= beginPos && EMaxIntvl[id].startT <= timePos) {//case 3
				continue;
			}
			else {//case 1 and 2
				int eLvalue = -eL[beginPos][id];
				
				EMaxIntvl[id].startT = beginPos;
				EMaxIntvl[id].endT = beginPos + eLvalue - 1;
				if (eLvalue < intvLen) {
					continue;
				}
				check = EMaxIntvl[id].endT - timePos;
				if (minEndT >= check) {
					minEndT = check;
					minEndTEdge = id;
				}
				EMaxIntvl[id].value = eW[beginPos][id];
				selectedNum++;
				edgeSetsR[check].emplace_back(id);
			}
		}
		//for (auto edgeIter = motifEdge->begin();
		//	edgeIter != edgeIterEnd; ++edgeIter) {
		//	id = edgeIter->id;
		//minEndTime[motifPos] = minEndT + intvE;
		//}
		edgeToMotifId->insert_or_assign(minEndTEdge, motifPos);
		//edgeToMotifId[minEndTEdge] = motifPos;
		delete tempMotif;
	}
}


#pragma region key incremental algorithms
/*DFTM (row number<=T-k+1)*/
void TGraphUEL::findTMotifsDynamic(int k,
	vec(TMotif*)*& newResult, int oriEndT,
	i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber, int TFchoice) {
#pragma region intialization
	i2iHMap vertex2Pos;//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp;
	//edgeToMotifId = DBG_NEW int[nEdge]; //the minimum right endpoint time of the interval of the motif can be expanded, records for each edge
	edgeToMotifId = DBG_NEW i2iHMap(); //the minimum right endpoint time of the interval of the motif can be expanded, records for each edge
	//minEndTime = DBG_NEW int[maxMotifNum]; //the minimum right endpoint time of the interval of the motif can be expanded, records for each edge
	motifSaveInfo = DBG_NEW SaveCCInfo[maxMotifNum]; //the SaveCCInfo of the motif which contains the edge

	//saveMotif = DBG_NEW bool[deltaEndT];//record motifNumSum[T] = motifNum[>T]
	int lineNum = endT - oriEndT + 1;
	vec(CComponents*) tempComponents;//temporary component list CC
	vec(int)* edgeSetsR = DBG_NEW vec(int)[lineNum];// R edge sets
	DisjointSet* disjointSet;//disjoint set

	veciter(int) iterEnd, iterStart;//iterator for edgeSetsR
	veciter(CComponents*) listEnd;//iterator for tempComponents
	int selectedNum;//the number of edges in edgeSetsR

	vec(int) saveMotifPos;//record the position of components which need to be saved
	/*record the position of components which need to combine with other components*/
	vec(int) combineMotifPos;

	int vertexNum;//the number of vertexs
	int realMotifNum;//the number of generated connected component
	int lastTime = oriEndT - k + 1;//last start time 
	int edgeEndT;//record the currently considering end time of edges
	veciter(TMotif*) resultEnd;
	int pos;
	int Te, intvE;//temporary end time
#pragma endregion

	//traverse all possible intervals' start time
	for (int Ts = startT; Ts <= lastTime; Ts++) {//O(T)
		selectedNum = 0;
		Te = Ts + k - 1;
		intvE = oriEndT;
		if (TFchoice == 2) {
			pos = ((Ts - startT) << 1) + 1;
		}
		else {
			pos = resultPos(Ts, intvE, startT, endT, k);
		}
		computeRESForDFTM(Ts, intvE,
			edgeSetsR, selectedNum,
			fixLabel, isEdgeTypeFixed, k,
			newResult, pos);//O(E)
		newResult[pos].clear();

		//initalize for every row
		disjointSet = DBG_NEW DisjointSet
		((int)((selectedNum << 1) + 1));
		root2Comp.clear();
		vertex2Pos.clear();
		vertexNum = 0;
		realMotifNum = 0;
		edgeEndT = endT;

		//R2L
		int startPos = endT - oriEndT;
		for (int tempT = startPos; tempT >= 0; tempT--, edgeEndT--) {//O(Es)
			//if (tempT >= maxIntervalLength) continue;
			iterStart = edgeSetsR[tempT].begin();
			iterEnd = edgeSetsR[tempT].end();
			if (iterStart == iterEnd) continue;
			//generateMaxTM and generatedExpTM
			genMotifInOneIntvDynamic(iterStart, iterEnd,
				vertex2Pos, disjointSet, vertexNum,
				combineMotifPos, realMotifNum, root2Comp, tempComponents,
				saveMotifPos, Ts, edgeEndT,
				newResult, k, motifNumber, TFchoice);
			edgeSetsR[tempT].clear();
		}

		//testing
		Test::updateMemoryUse();
		//Test::showMemoryUse();

		//release
		delete disjointSet;
		listEnd = tempComponents.end();
		for (auto listIter = tempComponents.begin();
			listIter != listEnd; ++listIter) {
			if (*listIter != NULL)
				delete *listIter;
		}
		tempComponents.clear();
	}
	delete[] edgeSetsR;
	//delete[] minEndTime;
	delete edgeToMotifId;
	delete[] motifSaveInfo;
}
#pragma endregion

void TGraphUEL::runDFTM(int k, vec(TMotif*)*& newResult, int oriEndT, i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber, int TFchoice) {
	//row number<=T-k+1
	findTMotifsDynamic(k, newResult,
		oriEndT, fixLabel, isEdgeTypeFixed, motifNumber, TFchoice);
	//row number>T-k+1  similar to FTM
	findTMotifs(k, newResult,
		fixLabel, isEdgeTypeFixed, motifNumber,
		oriEndT - k + 2, endT, TFchoice);
}