#pragma once
#include"stdafx.h"

class Edge{
public:
	Edge():s(-1),t(-1),id(-1){}

	Edge(int _s,int _t):id(-1) {
		if (_s < _t) {//ensure s < t
			s = _s;
			t = _t;
		}
		else {
			t = _s;
			s = _t;
		}
	}

	Edge(const Edge& e)
		:s(e.s), t(e.t), id(e.id) {}

	bool operator<(const Edge& edge1) const{
		return s < edge1.s || (s == edge1.s && t < edge1.t);
	}

	bool operator==(const Edge& edge1) const{
		return (s == edge1.s && t == edge1.t)||
			(s==edge1.t&&t==edge1.s);
	}

	friend ostream& operator<<(ostream& output, const Edge& e) {
		//output << "edge{startid: " << e.s << " endid: " << e.t << " }";
		output << e.s << "," << e.t;
		return output;
	}

	int s, t;//two endpoints, ensure s < t
	int id;//id of an edge
};
