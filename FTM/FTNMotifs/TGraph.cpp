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
	vec(TMotif*)*& result, int k, long long& motifNumber, int TFchoice) {
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
		motifEndT, motifNumber, TFchoice);
	Test::gne += clock() - begin;
}

void TGraph::genMotifInOneIntvDynamic(veciter(int)& iterStart, veciter(int)& iterEnd,
	i2iHMap& vertex2Pos, DisjointSet*& disjointSet, int& vertexNum,
	vec(int)& combineCCPos, int& realMotifNum,
	i2iHMap& root2Comp, vec(CComponents*)& tempComponents,
	vec(int)& saveCCPos, int motifStartT, int motifEndT,
	vec(TMotif*)*& result, int k, long long& motifNumber, int TFchoice) {
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
		motifEndT, motifNumber, TFchoice);
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

#pragma endregion 



void TGraph::generateExpTM(vec(int)&saveCCPos,
	vec(CComponents*)& tempComponents,
	vec(TMotif*)*& result, int k, int motifStartT, int motifEndT,
	long long& motifNumber, int TFchoice) {
	int savePos, tempSavePos;
	if (TFchoice == 2) {
		if (motifEndT != endT) savePos = (motifStartT - startT) << 1;
		else savePos = ((motifStartT - startT) << 1) + 1;
	}
	else {
		savePos = resultPos(motifStartT, motifEndT, startT, endT, k);
	}
	
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
			if (TFchoice == 2) {
				if (saveInfoIter->saveEndTime == endT) {
					tempSavePos = ((motifStartT - startT) << 1) + 1;
					linkMotif = result[tempSavePos][saveInfoIter->savePos];
					saveCCInfoP = motif->linkToMotifs(*saveInfoIter, linkMotif);
					saveinfo.saveCCInfoPos = saveCCInfoP;
					linkMotif->tempLinkToMotifs(saveinfo);
				}
				else {
					tempSavePos = (motifStartT - startT) << 1;
					linkMotif = result[tempSavePos][saveInfoIter->savePos];
					saveCCInfoP = motif->linkToMotifs(*saveInfoIter, linkMotif);
				}
			}
			else {
				tempSavePos = resultPos(motifStartT, saveInfoIter->saveEndTime, startT, endT, k);
				linkMotif = result[tempSavePos][saveInfoIter->savePos];
				saveCCInfoP = motif->linkToMotifs(*saveInfoIter, linkMotif);
				if (saveInfoIter->saveEndTime == endT) {
					saveinfo.saveCCInfoPos = saveCCInfoP;
					linkMotif->tempLinkToMotifs(saveinfo);
				}
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
	long long& motifNumber, int TFchoice) {
	int savePos, tempSavePos;
	if (TFchoice == 2) {
		if (motifEndT != endT) savePos = (motifStartT - startT) << 1;
		else savePos = ((motifStartT - startT) << 1) + 1;
	}
	else {
		savePos = resultPos(motifStartT, motifEndT, startT, endT, k);
	}
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
					if(TFchoice == 2)
						resultP = (motifStartT - startT) << 1;
					else
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
			if (TFchoice == 2) {
				if (saveInfoIter->saveEndTime == endT) {
					tempSavePos = ((motifStartT - startT) << 1) + 1;
					linkMotif = result[tempSavePos][saveInfoIter->savePos];
					saveCCInfoP = motif->linkToMotifs(*saveInfoIter, linkMotif);
					saveinfo.saveCCInfoPos = saveCCInfoP;
					linkMotif->tempLinkToMotifs(saveinfo);
				}
				else {
					tempSavePos = (motifStartT - startT) << 1;
					linkMotif = result[tempSavePos][saveInfoIter->savePos];
					saveCCInfoP = motif->linkToMotifs(*saveInfoIter, linkMotif);
				}
			}
			else {
				tempSavePos = resultPos(motifStartT, saveInfoIter->saveEndTime, startT, endT, k);
				linkMotif = result[tempSavePos][saveInfoIter->savePos];
				saveCCInfoP = motif->linkToMotifs(*saveInfoIter, linkMotif);
				if (saveInfoIter->saveEndTime == endT) {
					saveinfo.saveCCInfoPos = saveCCInfoP;
					linkMotif->tempLinkToMotifs(saveinfo);
				}
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

void TGraph::printMotif(vec(TMotif*)*& res, TMotif* motif, int k, int TFchoice) {
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
						int resultP;
						if (TFchoice == 2) {
							if (iter->saveEndTime != endT) {
								resultP = (motifStartT - startT) << 1;
							}
							else resultP = ((motifStartT - startT) << 1) + 1;
						}
						else {
							resultP = resultPos(motifStartT, iter->saveEndTime, startT, endT, k);
						}
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
