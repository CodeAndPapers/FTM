#include"TGraph.h"
#include"stdafx.h"

TGraph::TGraph(const TGraph& ances):
	startT(ances.startT), endT(ances.endT),
	nTimestamp(ances.nTimestamp),
	nNode(ances.nNode),nEdge(ances.nEdge),
	edge2ind(DBG_NEW map<Edge, int>(*ances.edge2ind)){
	edgeList = DBG_NEW Edge[nEdge];
	vertexToId = ances.vertexToId;
	ufset = DBG_NEW DisjointSet(*ances.ufset);
	for (int i = 0; i < nEdge; i++) {
		edgeList[i] = ances.edgeList[i];
	}
}

#pragma region save all information from temporal graph file
/*save all information from input temporal graph*/
void TGraph::loadInformation(const char*& src, set<Edge>& origEdge, 
	vec(int)& u_arr,vec(int)& v_arr,
	vec(int)& t_arr,vec(Label)& w_arr) {
	u_arr.reserve(ALLOC_MEM);
	v_arr.reserve(ALLOC_MEM);
	t_arr.reserve(ALLOC_MEM);
	w_arr.reserve(ALLOC_MEM);//allocate the memory

	nNode = nEdge = 0;
	//unordered_map<int, bool> nodeHashmap;//hashmap for recording node id 
	unordered_map<int, int>::iterator iter;
	int nodeNum = 0;//node num
	int firstT = 0x7fffffff, lastT = -1;//set startT and endT
	int flagT = -1;//check the same time as the first input data
	int u, v, t;
	Label w;
	edge2ind = DBG_NEW map<Edge, int>();// map from an edge to its index in edgeT

	FILE* file;
	file = fopen(src, "r+");
	if (!file) exit(0);
	char line[LINE_LENGTH];
	CLEARALL(line, 0, LINE_LENGTH, char);
	unordered_map<int, int> tempVertexToId;
	int sep1, sep2, sep3;//separator pos
	long ind = 0;
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
			if (firstT > t)firstT = t;
			if (lastT < t)lastT = t;
			if (flagT == -1) flagT = t;
			if (flagT == t) {//save static graph structure
				iter = tempVertexToId.find(u);
				if (iter == tempVertexToId.end()) {
					tempVertexToId[u] = nodeNum++;
				}//save the node id into hashmap
				iter = tempVertexToId.find(v);
				if (iter == tempVertexToId.end()) {
					tempVertexToId[v] = nodeNum++;
				}//save the map of name to id and id to name
				Edge e(u, v);
				//edgeSet->insert(e1);
				set<Edge>::iterator iter = origEdge.find(e);
				if (iter == origEdge.end()) {
					(*edge2ind)[e] = ind;
					e.id = ind;
					++ind;
					origEdge.insert(e);//save the edge of graph
				}
				else continue;//multiple edges
			}
			u_arr.emplace_back(u);
			v_arr.emplace_back(v);
			t_arr.emplace_back(t);
			w_arr.emplace_back(w);
		}
		fclose(file);
		this->startT = firstT;
		this->endT = lastT;
		this->nNode = nodeNum;
		this->ufset = DBG_NEW DisjointSet(nodeNum);
		this->nEdge = origEdge.size();
		this->nTimestamp = lastT - firstT + 1;
		this->edgeList = DBG_NEW Edge[nEdge];
		for (set<Edge>::iterator iter = origEdge.begin();
			iter != origEdge.end(); ++iter) {
			edgeList[iter->id] = *iter;
		}
	}
	catch (exception& e) {
		Util::printError(e);
		cout << "load error" << endl;
	}
}

/*save all information from input temporal graph with interval [fixedS, fixedE]*/
void TGraph::loadInformation(const char*& src, set<Edge>& origEdge,
	vec(int)& u_arr, vec(int)& v_arr,
	vec(int)& t_arr, vec(Label)& w_arr,
	int fixedS, int fixedE) {
	
	u_arr.reserve(ALLOC_MEM);
	v_arr.reserve(ALLOC_MEM);
	t_arr.reserve(ALLOC_MEM);
	w_arr.reserve(ALLOC_MEM);//allocate the memory

	nNode = nEdge = 0;
	unordered_map<int, int>::iterator iter;
	int nodeNum = 0;//node num
	int flagT = -1;//check the same time as the first input data
	int u, v, t;
	Label w;

	edge2ind = DBG_NEW map<Edge, int>();// map from an edge to its index in edgeT

	FILE* file;
	file = fopen(src, "r+");
	if (!file) exit(0);
	char line[LINE_LENGTH];
	CLEARALL(line,0, LINE_LENGTH, char);
	int sep1, sep2, sep3;//separator pos
	unordered_map<int, int> tempVertexToId;
	long ind = 0;
	if (fixedE == -1) fixedE = 0x7fffffff;
	int tempEndT = -1;
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
			if (flagT == -1) flagT = t;
			if (flagT == t) {//save static graph structure
				iter = tempVertexToId.find(u);
				if (iter == tempVertexToId.end()) {
					tempVertexToId[u] = nodeNum++;
				}//save the node id into hashmap
				iter = tempVertexToId.find(v);
				if (iter == tempVertexToId.end()) {
					tempVertexToId[v] = nodeNum++;
				}//save the map of name to id and id to name
				Edge e(u, v);
				set<Edge>::iterator iter = origEdge.find(e);
				if (iter == origEdge.end()) {
					(*edge2ind)[e] = ind;
					e.id = ind;
					++ind;
					origEdge.insert(e);//save the edge of graph
				}
				else continue;//multiple edges
			}
			if (t >= fixedS&&t <= fixedE) {
				if (t > tempEndT)tempEndT = t;
				u_arr.emplace_back(u);
				v_arr.emplace_back(v);
				t_arr.emplace_back(t);
				w_arr.emplace_back(w);
			}
		}
			
		fclose(file);
		this->startT = fixedS;
		this->endT = tempEndT;
		this->nNode = nodeNum;
		this->ufset = DBG_NEW DisjointSet(nodeNum);
		this->nEdge = origEdge.size();
		this->nTimestamp = fixedE - fixedS + 1;
		this->edgeList = DBG_NEW Edge[nEdge];
		set<Edge>::iterator iterEnd = origEdge.end();
		for (set<Edge>::iterator iter = origEdge.begin();
			iter != iterEnd; ++iter) {
			edgeList[iter->id] = *iter;
		}
	}
	catch (exception& e) {
		Util::printError(e);
		cout << "load error" << endl;
	}
}
#pragma endregion


#pragma region process of FTM

#pragma region FTM step 2 and step 3
/*generate motifs in one interval for FTM*/
void TGraph::genMotifInOneIntv(veciter(int)& iterStart, veciter(int)& iterEnd,
	i2iHMap& vertex2Pos, DisjointSet*& disjointSet, int& vertexNum,
	vec(int)& combineCCPos, int& realMotifNum,
	i2iHMap& root2Comp, vec(CComponents*)& tempComponents,
	vec(int)& saveCCPos, int motifStartT, int motifEndT,
	vec(TMotif*)*& result, int k, long long& motifNumber) {
	clock_t begin;
	begin = clock();

	#pragma region generateMaxTM
	/*fetch new edges from one R edge set and insert into the disjoint set (maintain connectivity)*/
	maintainUFSet(iterStart, iterEnd,
		vertex2Pos, disjointSet, vertexNum,
		combineCCPos);
	/*combine components which are connected after adding new edges from one R edge set,
		combine components before adding new edges for less computation 
		(combine components in generateMaxTM case 3)*/
	combineComponents(tempComponents,
		vertex2Pos, disjointSet, root2Comp,
		/*tempComponentsSize,*/ realMotifNum,
		combineCCPos);
	/*add edge into generated motif or generate new motif and check left expandable
		(generateMaxTM case 1, 2 and add edge into generated motif in generateMaxTM case 3)*/
	updateNewEdgeInfo(iterStart, iterEnd,
		tempComponents, vertex2Pos, disjointSet,
		root2Comp, /*tempComponentsSize,*/ realMotifNum, saveCCPos, motifStartT);//O(¦¤Es)
	Test::gm += clock()-begin;
	#pragma endregion

	begin = clock();
	generateExpTM(saveCCPos, tempComponents,
		result, k, motifStartT,
		motifEndT, motifNumber);
	Test::gne += clock() - begin;
}

void TGraph::genMotifInOneIntvDynamic(veciter(int)& iterStart, veciter(int)& iterEnd,
	i2iHMap& vertex2Pos, DisjointSet*& disjointSet, int& vertexNum,
	vec(int)& combineCCPos, int& realMotifNum,
	i2iHMap& root2Comp, vec(CComponents*)& tempComponents,
	vec(int)& saveCCPos, int motifStartT, int motifEndT,
	vec(TMotif*)*& result, int k, long long& motifNumber) {
	clock_t begin;
	begin = clock();

#pragma region generateMaxTM
	/*fetch new edges from one R edge set and insert into the disjoint set (maintain connectivity)*/
	maintainUFSet(iterStart, iterEnd,
		vertex2Pos, disjointSet, vertexNum,
		combineCCPos);
	/*combine components which are connected after adding new edges from one R edge set,
		combine components before adding new edges for less computation
		(combine components in generateMaxTM case 3)*/
	combineComponents(tempComponents,
		vertex2Pos, disjointSet, root2Comp,
		/*tempComponentsSize,*/ realMotifNum,
		combineCCPos);
	/*add edge into generated motif or generate new motif and check left expandable
		(generateMaxTM case 1, 2 and add edge into generated motif in generateMaxTM case 3)*/
	updateNewEdgeInfo(iterStart, iterEnd,
		tempComponents, vertex2Pos, disjointSet,
		root2Comp, /*tempComponentsSize,*/ realMotifNum, saveCCPos, motifStartT);//O(¦¤Es)
	Test::gm += clock() - begin;
#pragma endregion

	begin = clock();
	generateExpTMDynamic(saveCCPos, tempComponents,
		result, k, motifStartT,
		motifEndT, motifNumber);
	Test::gne += clock() - begin;
}


#pragma region generateMaxTM
/*fetch the new edge from one R edge set and insert into the disjoint set (maintain connectivity)*/
void TGraph::maintainUFSet(veciter(int)& infoBegin,
	veciter(int)& infoEnd, i2iHMap& vertex2Pos, DisjointSet*& ufset, int&vertexNum,
	vec(int)& combineCCPos) {
	//use disjoint set to maintain connectivity
	Edge* temp;
	int motifPos;
	i2iHMap_Iter vertexIt;
	int sId, tId, vertex;//vertexs'id in the disjoint set
	for (veciter(int) iter = infoBegin;
		iter != infoEnd; ++iter) {
		temp = &edgeList[*iter];//new edge from one R edge set
		vertex = temp->s;
		vertexIt = vertex2Pos.find(vertex);//O(1)
		if (vertexIt == vertex2Pos.end()) {
			sId = vertexNum;
			vertex2Pos[vertex] = vertexNum++;
		}
		else sId = vertexIt->second;
		vertex = temp->t;
		vertexIt = vertex2Pos.find(vertex);//O(1)
		if (vertexIt == vertex2Pos.end()) {
			tId = vertexNum;
			vertex2Pos[vertex] = vertexNum++;
		}
		else tId = vertexIt->second;

		/*union-find operation means that sId and tId are connected*/
		motifPos = ufset->newUnionVertexs(sId, tId);
		if (motifPos != -1)
			combineCCPos.push_back(motifPos);
	}
}

/*combine components which are connected after adding new edges from one R edge set,
		combine components before adding new edges for less computation
		(combine components in generateMaxTM case 3)*/
void TGraph::combineComponents(vec(CComponents*)& tempComponents,
	i2iHMap& vertex2Pos, DisjointSet*& disjointSet,
	i2iHMap& root2Comp, /*int& tempComponentsSize,*/int& realMotifNum,
	vec(int)& combineCCPos) {
#pragma region initialize
	int oldRoot, newRoot;//the original/current root in disjoint set
	int nowPos = 0;//component's new position in tempComponents
	CComponents* currentCComponents, *newCComponent;//move currentCComponents's edges to newCComponent
	
	vec(TEdge)* tempEdges;//edges of currentCComponents
	veciter(TEdge) edgeEnd;//tempEdges' iterator

	veciter(int) listEnd = combineCCPos.end();//tempComponents's iterator
	i2iHMap_Iter combineIter, mapIter;//root2Comp's iterator
	int combinePos;//the position of currentCComponents
#pragma endregion

#pragma region combination
	for (auto listIter = combineCCPos.begin();
		listIter != listEnd; ++listIter) {
		combineIter = root2Comp.find(*listIter);//cc position
		if (combineIter != root2Comp.end()) {//cc is combined to another cc
			combinePos = combineIter->second;
			currentCComponents = tempComponents[combinePos];
			
			oldRoot = currentCComponents->root;
			newRoot = disjointSet->find(oldRoot);
			root2Comp.erase(oldRoot);
			
			mapIter = root2Comp.find(newRoot);
			if (mapIter == root2Comp.end()) {//new edge will be added
				currentCComponents->root = newRoot;
				root2Comp[newRoot] = combinePos;
				continue;
			}

			tempEdges = &currentCComponents->edges;
			newCComponent = tempComponents[mapIter->second];
			edgeEnd = tempEdges->end();
			
			//insert edges of currentCComponents into newCComponent 
			for (auto edgeIter = tempEdges->begin();
				edgeIter != edgeEnd; ++edgeIter) {
				newCComponent->edges.emplace_back(*edgeIter);
			}
			vec(SaveCCInfo)* tempSaveInfo = &currentCComponents->saveInfo;//O(new insert edge + motif) <= O(original edges)
			veciter(SaveCCInfo) saveInfoEnd = tempSaveInfo->end();
			//insert saveInfo of currentCComponents into newCComponent 
			for (auto saveInfoIter = tempSaveInfo->begin();
				saveInfoIter != saveInfoEnd; ++saveInfoIter) {
				newCComponent->saveInfo.emplace_back(*saveInfoIter);
			}

			//maintain cc.intvl
			if (newCComponent->startT < currentCComponents->startT) {
				newCComponent->startT = currentCComponents->startT;
			}
			
			delete currentCComponents;//currentCComponents need to be deleted
			tempComponents[combinePos] = NULL;//this position will be deleted
			realMotifNum--;
		}
	}
#pragma endregion
	combineCCPos.clear();
}

/*add edge into generated motif or generate new motif and check left expandable
		(generateMaxTM case 1, 2 and add edge into generated motif in generateMaxTM case 3)*/
void TGraph::updateNewEdgeInfo(
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
	int startT;//the edge e's start time of EMaxIntvl[e]
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
		startT = EMaxIntvl[id].startT;//start time of the edge
		//startT = infoIter->startT;//start time of the edge
		label = EMaxIntvl[id].value; //label of the edge

		if (mapIter == root2Comp.end()) {//generateMaxTM case 1
			newCC = DBG_NEW CComponents(startT, root);
			//newCC->edges = DBG_NEW vector<TEdge>();
			newCC->edges.emplace_back(id, label);
			tempComponentsSize = (int)tempComponents.size();
			tempComponents.emplace_back(newCC);
			root2Comp[root] = tempComponentsSize;//update root2Comp

			if (startT == startTime) {//check left unexpandable
				saveCCPos.emplace_back(tempComponentsSize);
				hasSaved[tempComponentsSize] = true;
			}

			realMotifNum++;
		}
		else {//generateMaxTM case 2
			ccPos = mapIter->second;
			generatedCC = tempComponents[ccPos];
			generatedCC->edges.emplace_back(id, label);
			if (generatedCC->startT < startT) {
				generatedCC->startT = startT;
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
#pragma endregion 



#pragma region key incremental algorithms
/*DFTM (row number<=T-k+1)*/
void TGraph::findTMotifsDynamic(int k,
	vec(TMotif*)*& newResult, int oriEndT,
	i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber) {
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
		pos = resultPos(Ts, intvE, startT, endT, k);
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
				newResult, k, motifNumber);
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


#pragma region key static algorithms
/*FTM*/
void TGraph::findTMotifs(int k, vec(TMotif*)*& result,
	i2bHMap& fixLabel, bool isEdgeTypeFixed, long long& motifNumber,
	int choiceStartT, int choiceEndT) {
#pragma region initialize
	int lastTime = choiceEndT - k + 1;//last start time 
	i2iHMap vertex2Pos;//map the vertex's id to the position in disjoint set
	/*map the root in disjoint set to the position of its corresponding
	connected component in the component list*/
	i2iHMap root2Comp;

	vec(CComponents*) tempComponents;//temporary component list CC
	//int lineNum = lastTime - choiceStartT + 1;
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

void TGraph::generateExpTM(vec(int)&saveCCPos,
	vec(CComponents*)& tempComponents,
	vec(TMotif*)*& result, int k, int motifStartT, int motifEndT,
	long long& motifNumber) {
	int savePos = resultPos(motifStartT, motifEndT, startT, endT, k), tempSavePos;
	veciter(int) saveMotifEnd = saveCCPos.end();
	CComponents* tempCC;
	veciter(SaveCCInfo) saveCCInfoEnd;
	vec(SaveCCInfo)* saveCCInfoPtr;
	TMotif* motif,*linkMotif;
	int saveCCInfoP;
	//newly generated connected components
	/*if (motifStartT == 0 && motifEndT == 10) { 
		cout << endl;
	}*/
	for (auto saveMotifIter = saveCCPos.begin();
		saveMotifIter != saveMotifEnd; ++saveMotifIter) {
		tempCC = tempComponents[*saveMotifIter];
		motif = DBG_NEW TMotif(motifStartT, motifEndT);
		motif->setMotifEdge(tempCC->edges);
		 
		//update pointer to other motif(reuse other motif's edges)   O(motif number) <= O(original saved edge number)
		saveCCInfoPtr = &tempCC->saveInfo;
		saveCCInfoEnd = saveCCInfoPtr->end();
		SaveCCInfo saveinfo((int)result[savePos].size(), motifEndT);
		for (auto saveInfoIter = saveCCInfoPtr->begin();
			saveInfoIter != saveCCInfoEnd; ++saveInfoIter) {
			//reuse edges from original motifs
			tempSavePos = resultPos(motifStartT, saveInfoIter->saveEndTime, startT, endT, k);
			linkMotif = result[tempSavePos][saveInfoIter->savePos];
			saveCCInfoP = motif->linkToMotifs(*saveInfoIter, linkMotif);
			if (saveInfoIter->saveEndTime == endT) {
				saveinfo.saveCCInfoPos = saveCCInfoP;
				linkMotif->tempLinkToMotifs(saveinfo);
			}
		}

		tempCC->edges.clear();
		tempCC->saveInfo.clear();
		tempCC->saveInfo.emplace_back(saveinfo);

		result[savePos].emplace_back(motif);
	}
	motifNumber += saveCCPos.size();
	saveCCPos.clear();
}

void TGraph::generateExpTMDynamic(vec(int)&saveCCPos,
	vec(CComponents*)& tempComponents,
	vec(TMotif*)*& result, int k, int motifStartT, int motifEndT,
	long long& motifNumber) {
	int savePos = resultPos(motifStartT, motifEndT, startT, endT, k), tempSavePos;
	veciter(int) saveMotifEnd = saveCCPos.end();
	CComponents* tempCC;
	veciter(SaveCCInfo) saveCCInfoEnd;
	vec(SaveCCInfo)* saveCCInfoPtr;
	TMotif* motif, *linkMotif;
	int saveCCInfoP, resultP;
	int checkMotifId;
	i2iHMap_Iter mapEnd = edgeToMotifId->end(), mapIter;
	int edgeId;
	//newly generated connected components
	/*if (motifStartT == 0 && motifEndT == 10) {
		cout << endl;
	}*/
	for (auto saveMotifIter = saveCCPos.begin();
		saveMotifIter != saveMotifEnd; ++saveMotifIter) {
		tempCC = tempComponents[*saveMotifIter];
		motif = DBG_NEW TMotif(motifStartT, motifEndT);
		motif->setMotifEdge(tempCC->edges);

		//update pointer to this motif(reuse this motif's edges)   O(1)
		edgeId = tempCC->edges.rbegin()->id;
		mapIter = edgeToMotifId->find(edgeId);
		if (mapIter != mapEnd) {
			checkMotifId = mapIter->second;
			if (checkMotifId != -1) {
				SaveCCInfo& saveCCInfo = motifSaveInfo[checkMotifId];
				if (saveCCInfo.saveEndTime != -1) {//need to update the pointer
					resultP = resultPos(motifStartT, saveCCInfo.saveEndTime, startT, endT, k);
					linkMotif = result[resultP][saveCCInfo.savePos];
					saveCCInfoPtr = linkMotif->getOtherEdge();
					(*saveCCInfoPtr)[saveCCInfo.saveCCInfoPos] = SaveCCInfo(result[savePos].size(), motifEndT);
				}
				edgeToMotifId->insert_or_assign(edgeId, -1);
			}
		}

		//update pointer to other motif(reuse other motif's edges)   O(motif number) <= O(original saved edge number)
		saveCCInfoPtr = &tempCC->saveInfo;
		saveCCInfoEnd = saveCCInfoPtr->end();
		SaveCCInfo saveinfo((int)result[savePos].size(), motifEndT);
		for (auto saveInfoIter = saveCCInfoPtr->begin();
			saveInfoIter != saveCCInfoEnd; ++saveInfoIter) {
			//reuse edges from original motifs
			tempSavePos = resultPos(motifStartT, saveInfoIter->saveEndTime, startT, endT, k);
			linkMotif = result[tempSavePos][saveInfoIter->savePos];
			saveCCInfoP = motif->linkToMotifs(*saveInfoIter, linkMotif);
			if (saveInfoIter->saveEndTime == endT) {
				saveinfo.saveCCInfoPos = saveCCInfoP;
				linkMotif->tempLinkToMotifs(saveinfo);
			}
		}

		tempCC->edges.clear();
		tempCC->saveInfo.clear();
		tempCC->saveInfo.emplace_back(saveinfo);

		result[savePos].emplace_back(motif);
	}
	motifNumber += saveCCPos.size();
	saveCCPos.clear();
}
#pragma endregion

#pragma endregion

/*print motif*/
//void TGraph::printMotif(TMotif* motif) {
//	vec(TEdge)* motifEdge = motif->getMotifEdge();
//	veciter(TEdge) listEnd = motifEdge->end();
//	Edge *edge;
//	for (auto iter = motifEdge->begin();
//		iter != listEnd; ++iter) {
//		edge = &edgeList[iter->id];
//		cout << edge->s << "," << edge->t <<
//			"," << iter->label <<endl;
//	}
//}

void TGraph::printMotif(vec(TMotif*)*& res, TMotif* motif, int k) {
	int motifEndT = motif->getEndT(), motifStartT = motif->getStartT();
	//int counter = 0;
	if (motifEndT != endT) {
		queue<TMotif*> queue;
		queue.push(motif);
		TMotif* tempMotif;

		while (!queue.empty()) {
			tempMotif = queue.front();
			queue.pop();

			vec(TEdge)* motifEdge = tempMotif->getMotifEdge();
			veciter(TEdge) listEnd = motifEdge->end();
			Edge *edge;
			for (auto iter = motifEdge->begin();
				iter != listEnd; ++iter) {
				edge = &edgeList[iter->id];
				//counter++;
				cout << edge->s << "," << edge->t <<
					"," << iter->label << endl;
			}

			if (tempMotif->getEndT() != endT) {
				//reuse edges in other motifs
				vec(SaveCCInfo)* otherEdge = tempMotif->getOtherEdge();
				if (otherEdge != nullptr) {
					veciter(SaveCCInfo) otherListEnd = otherEdge->end();
					for (auto iter = otherEdge->begin();
						iter != otherListEnd; ++iter) {
						int resultP = resultPos(motifStartT, iter->saveEndTime, startT, endT, k);
						queue.push(res[resultP][iter->savePos]);
					}
				}
			}
		}
		
		/*if (motif->getSize() != counter) {
			cout << "!!!!!!!!"<< motifStartT << " " << motifEndT << " " << motif->getSize() << " " << counter << endl;
			exit(-1);
		}*/
	}
	else {
		vec(TEdge)* motifEdge = motif->getMotifEdge();
		veciter(TEdge) listEnd = motifEdge->end();
		Edge *edge;
		for (auto iter = motifEdge->begin();
			iter != listEnd; ++iter) {
			edge = &edgeList[iter->id];
			cout << edge->s << "," << edge->t <<
				"," << iter->label << endl;
		}
	}
}