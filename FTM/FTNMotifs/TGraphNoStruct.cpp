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
	CLEARALL(cE, -1, nEdge, int);

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
	cE = DBG_NEW int[nEdge];
	CLEARALL(cE, -1, nEdge, int);
	
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
	int choiceStartT, int choiceEndT, int TFchoice) {
#pragma region initialize
	int lastTime = choiceEndT - k + 1;//last start time 
	i2iHMap vertex2Pos;//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp;

	vec(CComponents*) tempComponents;//temporary component list CC
	vec(int)* edgeSetsR = DBG_NEW vec(int)[choiceEndT - choiceStartT + 1];// R edge sets
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

	for (int Ts = choiceEndT; Ts >= choiceStartT; Ts--) {//O(T)
		// select edges which keep fixed in [Ts,Ts+k-1]
		begin = clock();
		selectedNum = 0;
		computeRES(Ts, k,
			edgeSetsR, selectedNum,
			fixLabel, isEdgeTypeFixed);
		if (Ts > lastTime) { 
			Test::compr += clock() - begin;
			continue; 
		}

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
		Te = Ts + k - 1;
		//if (maxIntervalLength <= choiceEndT - Te) smt = 0;
		smt = edgeSetsR[choiceEndT - Te].size();
		maxSmt = max(maxSmt, smt);
		allSmt += smt;

		//R2L
		for (int tempT = choiceEndT - Te; tempT >= 0; tempT--, edgeEndT--) {
			//if (tempT >= maxIntervalLength) continue;
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
void TGraphNoStruct::updateNewEdgeInfo(
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
	int checkStartT;//check whether eW[currentP][id] == eW[currentP - 1][id]
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

		int currentP = startTime - startT;
		label = eW[currentP][id];//label of the edge
		if (currentP == 0 || label != eW[currentP - 1][id]) {
			checkStartT = startTime;
		}
		else {
			checkStartT = startTime - 1;
		}

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


void TGraphNoStruct::computeRES(int scanT, int k,
	vec(int)*& edgeSetsR, int& selectedNum,
	unordered_map<Label, bool>& fixLabel,
	bool isEdgeTypeFixed) {

	int scanPos = scanT - startT;
	Label edgeType;
	unordered_map<int, bool>::iterator fixLabelEnd = fixLabel.end();
	for (int i = 0; i < nEdge; i++) {
		edgeType = eW[scanPos][i];
		/*if (isEdgeTypeFixed){
			if(fixLabel.find(edgeType) == fixLabelEnd)
				continue;
		}*/

		if (cE[i] == -1 || edgeType != eW[cE[i]][i]) {
			cE[i] = scanPos;
		}
		int check = cE[i] - (scanPos + k - 1);
		if (check >= 0) {
			selectedNum++;
			edgeSetsR[check].emplace_back(i);
		}
	}
}

/*compRES for DFTM (row number<=T-k+1)*/
void TGraphNoStruct::computeRESForDFTM(int scanT, int oriEndT, int k,
	vec(int)*& edgeSetsR, int& selectedNum,
	unordered_map<Label, bool>& fixLabel,
	bool isEdgeTypeFixed, vec(TMotif*)*& allResult, int pos) {
	vec(TMotif*)& result = allResult[pos];
	int scanPos = scanT - startT;
	int timePos = oriEndT - startT;
	int check;
	Label edgeType;
	unordered_map<int, bool>::iterator fixLabelEnd = fixLabel.end();
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
		iter != iterEnd; ++iter, ++motifPos) {//consider edges only in TF^o[m,T]
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

			edgeType = eW[scanPos][id];
			/*if (isEdgeTypeFixed){
				if(fixLabel.find(edgeType) == fixLabelEnd)
					continue;
			}*/

			check = cE[id] - timePos;
			if (check >= 0) {
				if (minEndT >= check) {
					minEndT = check;
					minEndTEdge = id;
				}
				selectedNum++;
				edgeSetsR[check].emplace_back(id);
			}
		}
		
		edgeToMotifId->insert_or_assign(minEndTEdge, motifPos);
		delete tempMotif;
	}
}

/*get the label of edge with edgeId at the time*/
Label TGraphNoStruct::getWeight(int time, int edgeId) {
	return eW[time - startT][edgeId];
}


#pragma region key incremental algorithms
/*DFTM (row number<=T-k+1)*/
void TGraphNoStruct::findTMotifsDynamic(int k,
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
	int lastTime = oriEndT - k + 1, realLastTime = endT - k + 1;//last start time 
	int edgeEndT;//record the currently considering end time of edges
	veciter(TMotif*) resultEnd;
	int pos;
	int intvE;//temporary end time
#pragma endregion

	//traverse all possible intervals' start time
	for (int Ts = lastTime; Ts >= startT; Ts--) {//O(T)
		selectedNum = 0;
		intvE = oriEndT;
		if (TFchoice == 2) {
			pos = ((Ts - startT) << 1) + 1;
		}
		else {
			pos = resultPos(Ts, intvE, startT, endT, k);
		}
		computeRESForDFTM(Ts, oriEndT, k,
			edgeSetsR, selectedNum,
			fixLabel, isEdgeTypeFixed, 
			newResult, pos);//O(E)
		newResult[pos].clear();
		if (Ts > realLastTime) {
			continue;
		}
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


void TGraphNoStruct::runDFTM(int k, vec(TMotif*)*& newResult, int oriEndT, i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber, int TFchoice) {
	//row number>T-k+1  similar to FTM
	findTMotifs(k, newResult,
		fixLabel, isEdgeTypeFixed, motifNumber,
		oriEndT - k + 2, endT, TFchoice);
	//row number<=T-k+1
	findTMotifsDynamic(k, newResult,
		oriEndT, fixLabel, isEdgeTypeFixed, motifNumber, TFchoice);
}