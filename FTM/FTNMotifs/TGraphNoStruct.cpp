#include "TGraphNoStruct.h"
#include "stdafx.h"

TGraphNoStruct::TGraphNoStruct(const TGraphNoStruct& ances) :TGraph(ances) {
	eW = DBG_NEW Label*[nTimestamp];
	for (int i = 0; i < nTimestamp; i++) {
		eW[i] = DBG_NEW Label[nEdge];
		for (int j = 0; j < nEdge; j++) {
			eW[i][j] = ances.eW[i][j];
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
void TGraphNoStruct::constructGraph(const char* src) {
	vec(int) u_arr, v_arr, t_arr;
	vec(Label) w_arr;
	set<Edge> origEdge;
	loadInformation(src, origEdge, u_arr, v_arr, t_arr, w_arr);
	createStructure(u_arr, v_arr, t_arr, w_arr/*, name2id*/);
}

/*construct the temporal graph, but fixed the startT and endT*/
void TGraphNoStruct::constructGraph(const char* src,
	int startT, int endT) {
	vec(int) u_arr, v_arr, t_arr;
	vec(Label) w_arr;
	set<Edge> origEdge;
	loadInformation(src, origEdge, u_arr, v_arr, t_arr, w_arr,
		startT, endT);
	createStructure(u_arr, v_arr, t_arr, w_arr/*, name2id*/);

}

//increase snapshots of the temporal graph
void TGraphNoStruct::changeGraph(const char* src, int oldEndT, int limitNewEndT) {
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
	delete[] EMaxIntvl;
	EMaxIntvl = DBG_NEW Intv[nEdge];
	for (int j = 0; j < nEdge; j++) {
		EMaxIntvl[j].endT = EMaxIntvl[j].startT = -1;
	}
	for (int i = 0; i < nTimestamp; i++) {
		neweW[i] = DBG_NEW Label[nEdge];
		if (i < oldTimestamp) {
			for (int j = 0; j < nEdge; j++) {
				neweW[i][j] = eW[i][j];
			}
			delete[] eW[i];
		}
		else {
			for (int j = 0; j < nEdge; j++) {
				neweW[i][j] = 0x7fffffff;
			}
		}
	}
	delete[] eW;
	eW = neweW;
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
void TGraphNoStruct::createStructure(vec(int)& u_arr, vec(int)& v_arr,
	vec(int)& t_arr, vec(Label)& w_arr) {
	eW = DBG_NEW Label*[nTimestamp];
	EMaxIntvl = DBG_NEW Intv[nEdge];
	for (int j = 0; j < nEdge; j++) {
		EMaxIntvl[j].endT = EMaxIntvl[j].startT = -1;
	}
	for (int i = 0; i < nTimestamp; i++) {
		eW[i] = DBG_NEW Label[nEdge];
		for (int j = 0; j < nEdge; j++) {
			eW[i][j] = 0x7fffffff;
		}
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
		Edge e1(u, v);
		edgeInd = edge2ind->find(e1)->second;

		pos = t - startT;
		if (eW[pos][edgeInd] == 0x7fffffff) {
			eW[pos][edgeInd] = w;
		}
	}
	cout << "no struct: " << 4 * nTimestamp * nEdge << "Byte" << endl;
	endTime = clock();
	cout << "load: " << endTime - startTime << "ms" << endl;
}


#pragma region key static algorithms
/*FTM*/
void TGraphNoStruct::findTMotifs(int k, vec(TMotif*)*& result,
	i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber,
	int choiceStartT, int choiceEndT) {
//#pragma region initialize
//	int lastTime = choiceEndT - k + 1;//last start time 
//	i2iHMap vertex2Pos;//map the vertex's id to the position in disjoint set
//	/*map the root in disjoint set to the position of its corresponding
//	connected component in the component list*/
//	i2iHMap root2Comp;
//
//	vec(CComponents*) tempComponents;//temporary component list CC
//	int lineNum = lastTime - choiceStartT + 1;
//	//int setsRNum = (lineNum * (lineNum + 1)) >> 1;
//	//SAVEINFONOST_Vec* edgeSetsR = DBG_NEW SAVEINFONOST_Vec[setsRNum];// R edge sets
//	SAVEINFO_Vec* edgeSetsR = DBG_NEW SAVEINFO_Vec[lineNum];// R edge sets
//	DisjointSet* disjointSet;//disjoint set
//
//	//SAVEINFONOST_VIter iterEnd, iterStart;//iterator for edgeSetsR
//	SAVEINFO_VIter iterEnd, iterStart;//iterator for edgeSetsR
//	veciter(CComponents*) listEnd;//iterator for tempComponents
//	int selectedNum;//the number of edges in edgeSetsR
//
//	vec(int) saveMotifPos;//record the position of components which need to be saved
//	/*record the position of components which need to combine with other components*/
//	vec(int) combineMotifPos;
//	long long maxSelect = 0;//the max number of edge in G(S[m,m+k-1])
//	long long selectedSum = 0;//the sum of edges in G(S[m,m+k-1])
//	clock_t begin;//record time
//	long long allSmt = 0, maxSmt = 0, smt;//the sum and the max number of S[m,T]
//
//	int vertexNum;//the number of vertexs
//	int realMotifNum;//the number of generated connected component
//	int edgeEndT;//record the currently considering end time of R edge sets
//	int Te;//the end time of interval
//#pragma endregion
//	//int tempPosForEdgeSetsR = 0;
//	int *scanP = DBG_NEW int[nEdge];
//	CLEARALL(scanP, 0, nEdge, int);
//	for (int Ts = choiceStartT; Ts <= lastTime; Ts++/*, tempPosForEdgeSetsR += lineNum, lineNum--*/) {//O(T)
//		// select edges which keep fixed in [Ts,Ts+k-1]
//		begin = clock();
//		selectedNum = 0;
//		Te = Ts + k - 1;
//		computeRES(Ts, Te, scanP,
//			edgeSetsR, selectedNum,
//			fixLabel, isEdgeTypeFixed);
//
//		//initalize for every row
//		maxSelect = max(maxSelect, selectedNum);
//		selectedSum += selectedNum;
//		disjointSet = DBG_NEW DisjointSet((int)((selectedNum << 1) + 1));
//		vertex2Pos.clear();
//		root2Comp.clear();
//		vertexNum = 0;
//		realMotifNum = 0;
//		edgeEndT = choiceEndT;
//		Test::compr += clock() - begin;
//
//		//Test S(m,T)
//		//smt = edgeSetsR[tempPosForEdgeSetsR + choiceEndT - Te].size();
//		smt = edgeSetsR[choiceEndT - Te].size();
//		maxSmt = max(maxSmt, smt);
//		allSmt += smt;
//
//		//R2L
//		for (int tempT = choiceEndT - Te; tempT >= 0; tempT--, edgeEndT--) {
//			iterStart = edgeSetsR[tempT].begin();
//			iterEnd = edgeSetsR[tempT].end();
//			//iterStart = edgeSetsR[tempPosForEdgeSetsR + tempT].begin();
//			//iterEnd = edgeSetsR[tempPosForEdgeSetsR + tempT].end();
//			if (iterStart == iterEnd) continue;
//			
//			//generateMaxTM and generatedExpTM
//			genMotifInOneIntv(iterStart, iterEnd,
//				vertex2Pos, disjointSet, vertexNum,
//				combineMotifPos, realMotifNum, root2Comp, tempComponents,
//				saveMotifPos, Ts, edgeEndT,
//				result, k, motifNumber);
//
//			edgeSetsR[tempT].clear();
//		}
//
//		//testing
//		Test::updateMemoryUse();
//		//Test::showMemoryUse();
//
//		//release
//		delete disjointSet;
//		listEnd = tempComponents.end();
//		for (veciter(CComponents*) listIter = tempComponents.begin();
//			listIter != listEnd; ++listIter) {
//			if (*listIter != NULL)
//				delete *listIter;
//		}
//		tempComponents.clear();
//	}
//	delete[] edgeSetsR;
//	delete[] scanP;
//	cout << SELECT_EDGE << maxSelect << endl;
//	cout << MEAN_SELECT_EDGE << selectedSum / (lastTime - choiceStartT + 1) << endl;
//	cout << "maxSmT: " << maxSmt << " averSmT: " << allSmt / (lastTime - choiceStartT + 1) << endl;
}
#pragma endregion

/*generate motifs in one interval for FTM*/
//void TGraphNoStruct::genMotifInOneIntv(SAVEINFONOST_VIter& iterStart, SAVEINFONOST_VIter& iterEnd,
//	i2iHMap& vertex2Pos, DisjointSet*& disjointSet, int& vertexNum,
//	vec(int)& combineCCPos, int& realMotifNum,
//	i2iHMap& root2Comp, vec(CComponents*)& tempComponents,
//	vec(int)& saveCCPos, int motifStartT, int motifEndT,
//	vec(TMotif*)*& result, int k, long long& motifNumber) {
//	clock_t begin;
//	begin = clock();
//
//#pragma region generateMaxTM
//	/*fetch new edges from one R edge set and insert into the disjoint set (maintain connectivity)*/
//	maintainUFSet(iterStart, iterEnd,
//		vertex2Pos, disjointSet, vertexNum,
//		combineCCPos);
//	/*combine components which are connected after adding new edges from one R edge set,
//		combine components before adding new edges for less computation
//		(combine components in generateMaxTM case 3)*/
//	combineComponents(tempComponents,
//		vertex2Pos, disjointSet, root2Comp,
//		/*tempComponentsSize,*/ realMotifNum,
//		combineCCPos);
//	/*add edge into generated motif or generate new motif and check left expandable
//		(generateMaxTM case 1, 2 and add edge into generated motif in generateMaxTM case 3)*/
//	updateNewEdgeInfo(iterStart, iterEnd,
//		tempComponents, vertex2Pos, disjointSet,
//		root2Comp, /*tempComponentsSize,*/ realMotifNum, saveCCPos, motifStartT);//O(¦¤Es)
//	Test::gm += clock() - begin;
//#pragma endregion
//
//	begin = clock();
//	generateExpTM(saveCCPos, tempComponents,
//		result, k, motifStartT,
//		motifEndT, motifNumber);
//	Test::gne += clock() - begin;
//}
//
//
//#pragma region generateMaxTM
///*fetch the new edge from one R edge set and insert into the disjoint set (maintain connectivity)*/
//void TGraphNoStruct::maintainUFSet(SAVEINFONOST_VIter& infoBegin,
//	SAVEINFONOST_VIter& infoEnd, i2iHMap& vertex2Pos, DisjointSet*& ufset, int&vertexNum,
//	vec(int)& combineCCPos) {
//	//use disjoint set to maintain connectivity
//	Edge* temp;
//	int motifPos;
//	i2iHMap_Iter vertexIt;
//	int sId, tId, vertex;//vertexs'id in the disjoint set
//	for (SAVEINFONOST_VIter iter = infoBegin;
//		iter != infoEnd; ++iter) {
//		temp = &edgeList[iter->edgeId];//new edge from one R edge set
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
//		if (motifPos != -1)
//			combineCCPos.push_back(motifPos);
//	}
//}
//
///*add edge into generated motif or generate new motif and check left expandable
//		(generateMaxTM case 1, 2 and add edge into generated motif in generateMaxTM case 3)*/
//void TGraphNoStruct::updateNewEdgeInfo(
//	SAVEINFONOST_VIter& infoBegin, SAVEINFONOST_VIter& infoEnd,
//	vec(CComponents*)& tempComponents,
//	i2iHMap& vertex2Pos, DisjointSet*& disjointSet,
//	i2iHMap& root2Comp, /*int& tempComponentsSize,*/ int& realMotifNum,
//	vec(int)& saveCCPos, int startTime) {
//
//#pragma region initialize
//	i2bHMap hasSaved;//means whether the motif has been saved 
//	i2iHMap_Iter mapIter;//root2Comp's iterator 
//	int root;//the root of one vertex in the disjoint set
//	int id;//the edge e's id
//	int startT;//the edge e's start time of EMaxIntvl[e]
//	int label;//the edge e's label
//	int ccPos;//cc's position in tempComponents
//	/*the cc which has already generated/will generate*/
//	CComponents* generatedCC, *newCC;
//	int tempComponentsSize;//size of tempComponents
//#pragma endregion
//
//	for (SAVEINFONOST_VIter infoIter = infoBegin;
//		infoIter != infoEnd; ++infoIter) {//new edges in R edge sets
//
//		id = infoIter->edgeId;//edge's id
//
//		/*the root of the edge's vertex in the disjoint set*/
//		root = disjointSet->find(vertex2Pos[edgeList[id].s]);
//
//		mapIter = root2Comp.find(root);//check which cc the edge belongs to 
//		startT = infoIter->startT;//start time of the edge
//		label = infoIter->label; //label of the edge
//
//		if (mapIter == root2Comp.end()) {//generateMaxTM case 1
//			newCC = DBG_NEW CComponents(startT, root);
//			//newCC->edges = DBG_NEW vector<TEdge>();
//			newCC->edges.emplace_back(id, label);
//			tempComponentsSize = (int)tempComponents.size();
//			tempComponents.emplace_back(newCC);
//			root2Comp[root] = tempComponentsSize;//update root2Comp
//
//			if (startT == startTime) {//check left unexpandable
//				saveCCPos.emplace_back(tempComponentsSize);
//				hasSaved[tempComponentsSize] = true;
//			}
//
//			realMotifNum++;
//		}
//		else {//generateMaxTM case 2
//			ccPos = mapIter->second;
//			generatedCC = tempComponents[ccPos];
//			generatedCC->edges.emplace_back(id, label);
//			if (generatedCC->startT < startT) {
//				generatedCC->startT = startT;
//			}
//			/* check whether this cc has already been inserted into save list */
//			if (hasSaved.find(ccPos) == hasSaved.end()) {
//				if (generatedCC->startT == startTime) {//check left unexpandable
//					saveCCPos.emplace_back(ccPos);
//					hasSaved[ccPos] = true;
//				}
//			}
//		}
//	}
//}
//#pragma endregion 

void TGraphNoStruct::computeRES(int choiceStartT, int choiceEndT,
	SAVEINFONOST_Vec*& edgeSetsR,
	unordered_map<Label, bool>& fixLabel,
	bool isEdgeTypeFixed, int k, int setsRNum) {

	int lastTimePos = choiceEndT - startT;
	int firstTimePos = choiceStartT - startT;
	int edgeSetsRPosChange = choiceEndT - k + 1 - startT;
	Label preEdgeType, curEdgeType;
	unordered_map<int, bool>::iterator fixLabelEnd = fixLabel.end();
	int preSameLabelPos;
	int edgeSetsRPos, removePos;
	bool saveToRSet;
	int preSavePForOneSet, preSaveP;
	for (int i = 0; i < nEdge; i++) {
		int timePos = lastTimePos;
		preSameLabelPos = timePos;
		preEdgeType = eW[timePos][i];
		timePos--;
		saveToRSet = false;
		
		edgeSetsRPos = setsRNum;
		removePos = 1;
		
		for (; timePos >= firstTimePos; timePos--) {
			if (timePos <= edgeSetsRPosChange) {
				edgeSetsRPos -= removePos;
				removePos++;
			}
			curEdgeType = eW[timePos][i];
			if (preEdgeType != curEdgeType) {
				preSameLabelPos = timePos;
				preEdgeType = curEdgeType;
				saveToRSet = false;
			}
			else if(saveToRSet){
				edgeSetsR[preSaveP][preSavePForOneSet].startT--; //expandable
			}

			if (preSameLabelPos - timePos + 1 >= k && (!isEdgeTypeFixed || fixLabel.find(curEdgeType) != fixLabelEnd)) {
				saveToRSet = true;
				preSaveP = edgeSetsRPos + preSameLabelPos - (timePos + k - 1);
				edgeSetsR[preSaveP].emplace_back(i,
					timePos, //initialize
					curEdgeType);
				preSavePForOneSet = edgeSetsR[preSaveP].size() - 1;
			}
		}
	}
}

void TGraphNoStruct::computeRES(int intvB, int intvE, int*& scanP,
	SAVEINFO_Vec*& edgeSetsR, int& selectedNum,
	unordered_map<Label, bool>& fixLabel,
	bool isEdgeTypeFixed) {
	int timePos = intvE - startT;
	int beginPos = intvB - startT;
	int check;
	Label edgeType;
	unordered_map<int, bool>::iterator fixLabelEnd = fixLabel.end();
	int intvLen = intvE - intvB + 1;
	int lastPos = endT - startT;
	int intvStart;
	//if (maxeL[timePos] < intvLen) return;
	for (int i = 0; i < nEdge; i++) {
		edgeType = eW[beginPos][i];
		/*not exists the maximum interval containing [intvB,intvE] for case 1 and 2
			put here for less time*/
		if ((isEdgeTypeFixed && fixLabel.find(edgeType)
				== fixLabelEnd)) {
			continue;
		}

		intvStart = scanP[i];
		if (intvStart > timePos) {
			selectedNum++;
			check = EMaxIntvl[i].endT - timePos;
			edgeSetsR[check].emplace_back(i, edgeType);
		}
		else if (intvStart > beginPos) {
			continue;
		}
		else {//intvStart == beginPos
			int j = intvStart + 1;
			for (; j <= lastPos + 1; j++) {
				if (j == lastPos + 1) {
					selectedNum++;
					EMaxIntvl[i].startT = intvStart;
					EMaxIntvl[i].endT = lastPos;
					check = lastPos - timePos;
					edgeSetsR[check].emplace_back(i, edgeType);
					scanP[i] = j;
				}
				else if (eW[j][i] != edgeType) {
					if (j - 1 >= timePos) {
						selectedNum++;
						check = j - timePos - 1;
						EMaxIntvl[i].startT = intvStart;
						EMaxIntvl[i].endT = j - 1;
						edgeSetsR[check].emplace_back(i, edgeType);
					}
					scanP[i] = j;
					break;
				}
			}
		}
	}
}

/*get the label of edge with edgeId at the time*/
Label TGraphNoStruct::getWeight(int time, int edgeId) {
	return eW[time - startT][edgeId];
}

