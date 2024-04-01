#pragma once
#include "stdafx.h"
#define LINE_LENGTH 1000
#define SEP_CHAR ','
#define STR2INT(a) Util::stringToInt(a)
#define STR2BOOL(a) Util::stringToBool(a)

using namespace std;
class Util {
public:
	
	//transform char* to int
	inline static int stringToInt(char* str) {
		return atoi(str);
	}

	//transform char* to bool
	inline static bool stringToBool(char* str) {
		return 0!=atoi(str);
	}
	
	//print exception
	inline static void printError(const exception& e) {
		cout <<"exception caught:"<< e.what() << endl;
	}
	
	//return the maximum number among a,b,c
	inline static int getMax(int a,int b,int c) {
		int max = a > b ? a : b;
		return max > c ? max : c;
	}
	
	//return the minimum number among a,b,c
	inline static int getMin(int a, int b, int c) {
		int min = a < b ? a : b;
		return min < c ? min : c;
	}

};


#pragma region testing
class Test {
public:
	static long long sumIntvLen, maxIntvLen;//used for testing
	//static unordered_map<int, int> hist;
	//static unordered_map<int, int> vertexHist;
	static clock_t compr, gm, gne;

	static SIZE_T peakMemory;

	static void updateMemoryUse() {
		//get the handle of the current process
		HANDLE currentProcess = GetCurrentProcess();
		PROCESS_MEMORY_COUNTERS pmc;
		GetProcessMemoryInfo(currentProcess, &pmc, sizeof(pmc));
		if (peakMemory < pmc.WorkingSetSize) peakMemory = pmc.WorkingSetSize;
	}

	static void showRealPeakMemoryUse() {
		cout << "Real Peak Memory Use:" <<
			peakMemory << "Byte" << endl;
	}

	//print the working set size in bytes
	static void showMemoryUse() {
		//get the handle of the current process
		HANDLE currentProcess = GetCurrentProcess();
		PROCESS_MEMORY_COUNTERS pmc;
		GetProcessMemoryInfo(currentProcess, &pmc, sizeof(pmc));
		cout << "Memory Use:" <<
			pmc.WorkingSetSize << "Byte" << endl;
	}

	//print the peak working set size in bytes
	static void showPeakMemoryUse() {
		//get the handle of the current process
		HANDLE currentProcess = GetCurrentProcess();
		PROCESS_MEMORY_COUNTERS pmc;
		GetProcessMemoryInfo(currentProcess, &pmc, sizeof(pmc));
		cout << "Peak Memory Use:" <<
			pmc.PeakWorkingSetSize << "Byte" << endl;
	}
};
#define TIMES_PER_SEC (1.0e9)
#pragma endregion