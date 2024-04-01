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
		if (isEdgeTypeFixed){
			edgeType = eW[beginPos][i];
			if(fixLabel.find(edgeType) == fixLabelEnd) 
				continue;
		}
		
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
	int lastPos = endT - startT;
	int /*maxIter=0,*/intvStart;
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
		delete tempMotif;
	}
}
