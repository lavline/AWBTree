#ifndef _DATA_STRUCTURE_H
#define _DATA_STRUCTURE_H
#include <string>
#include <vector>
#include <stdint.h>

using namespace std;

struct IntervalPred{
    int att;
    int lowVal, highVal;
};

struct IntervalSub{
    int id;
    int size;
    vector<IntervalPred> pred;
};

struct Element{
    int lowVal, highVal;
    int subID;
    bool operator!=(const Element b) const
    {
        return this->subID != b.subID;
    }
};

struct Pair{
    int att;
    int value;
};

struct Pub{
    int size;
    vector<Pair> pairs;
};


#endif //_DATA_STRUCTURE_H
