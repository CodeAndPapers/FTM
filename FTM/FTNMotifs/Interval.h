#pragma once
#include "stdafx.h"
/* interval with start time, end time and additional information (label)*/

template <class Scalar, typename Info>
class Interval {
public:
	//guarantee the start<=end
	Interval() :startT(-1), endT(-1),value(-1) {}
	
	/*Interval(const Scalar& start, const Scalar& end, Info val)
		:startT(start),endT(end),value(val){}

	Interval(const Scalar& start, const Scalar& end)
		:startT(start), endT(end), value(-1) {}*/

	Interval(const Scalar& start, const Scalar& end, Info val)
		:startT(start), endT(end), value(val) {}

	Interval(const Interval & intv)
		:startT(intv.startT), endT(intv.endT), value(intv.value) {}

	void setValue(const Scalar& start, const Scalar& end, Info val) {
		this->startT = start;
		this->endT = end;
		this->value = val;
	}

	bool operator<(const Interval & intv) const {
		return (startT < intv.startT ||
			(startT == intv.startT && endT < intv.endT) );
	}

	bool operator>(const Interval & intv) const {
		return (startT > intv.startT ||
			(startT == intv.startT && endT > intv.endT) );
	}

	bool operator<=(const Interval & intv) const {
		return (startT <= intv.startT ||
			(startT == intv.startT && endT <= intv.endT) );
	}

	bool operator>=(const Interval & intv) const {
		return (startT >= intv.startT ||
			(startT == intv.startT && endT >= intv.endT) );
	}

	bool operator==(const Interval & intv) const {
		return (startT == intv.startT && endT == intv.endT
			&& value == intv.value
			&& edgeId=intv.edgeId);
	}

	friend ostream& operator<<(ostream& output, const Interval & e) {
		output << "Interval{" << e.startT <<
		   "," << e.endT << ":" << e.value << "}"<<endl;
		return output;
	}

	Scalar startT, endT;//start time and end time
	Info value;//additional information: label
};
