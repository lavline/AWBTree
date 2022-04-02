#include <iostream>
//#include <random>
//#include "chrono_time.h"
#include "AWBTree.h"

using namespace std;

int main()
{
    /*std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 10000);*/
    int totalsubs = 1000000;
    int subs = SUB_NUM;     // number of subscriptions.
    int pubs = 1000;     // number of publications.
    int atts = 40;     // total number of attributes, i.e. dimensions.
    int cons = 10;     // number of constraints(predicates) in one sub.
    int m = 40;        // number of constraints in one pub.
    int attDis = 0;   // the distribution of attributes in subs and pubs. 0:uniform distribution | 1:zipf distribution
    int valDis = 0;   // the distribution of values in subs and pubs. 0:uniform, 1: both zipf, -1: sub zipf but pub uniform
    int valDom = 1000000;   // cardinality of values.
    double alpha = 0; // parameter for zipf distribution.
    double width = 0.7; // width of a predicate.
    int genRand = 1;
    float dPoint = 0.2;//partition point
    uint16_t branch = 2000;
    uint16_t w_size = 40; //width cell size
    
    AWBTree a(atts, m, subs, pubs, w_size, branch, valDom, dPoint);

    string dir = "/data/hq/.vs/genData-1.0/";
    string sub_file = dir + "sub-" + to_string(totalsubs / 1000) + "K-" + to_string(atts) + "D-" + to_string(cons) + "A-" + to_string(attDis) + "Ad-"
        + to_string(valDis) + "Vd-" + to_string((int)(alpha * 10)) + "al-" + to_string((int)(width * 10)) + "W-" + to_string(genRand) + "R.txt";
    string pub_file = dir + "pub-" + to_string(atts) + "D-" + to_string(m) + "A-" + to_string(attDis) + "Ad-"
        + to_string(valDis) + "Vd-" + to_string((int)(alpha * 10)) + "al.txt";
    cout << sub_file << endl;
    cout << pub_file << endl;

    if (!a.readSubs(sub_file)) return -1;
    if (!a.readPus(pub_file)) return -1;

    vector<double> insertTimeList;
    vector<double> matchTimeList;
    vector<double> deleteTimeList;
    vector<uint32_t> matchSubList;

    //insert PobaTree
    cout << "insert start...\n";

    for (int i = 0; i < subs; i++)
    {
        Timer subStart;
        a.insert(a.subList[i]); // Insert sub[i] into data structure.
        int64_t insertTime = subStart.elapsed_nano(); // Record inserting time in nanosecond.
        insertTimeList.push_back((double)insertTime / 1000000);
    }
    a.setbits();
    double avgInsertTime = 0;
    for (int i = 0; i < insertTimeList.size(); i++)avgInsertTime += insertTimeList[i];
    avgInsertTime /= insertTimeList.size();

    cout << "insert complete! avgInsertTime = " << avgInsertTime << "ms\n";

    double mem = a.memory();
    cout << "memory cost " << mem / 1024 / 1024 << "MB\n";

    // match
    for (int i = 0; i < pubs; i++)
    {
        int matchSubs = 0;                              // Record the number of matched subscriptions.
        Timer matchStart;

        //a.forward(a.pubList[i], matchSubs);
        //a.backward(a.pubList[i], matchSubs);
        //a.hybrid(a.pubList[i], matchSubs);

        int64_t eventTime = matchStart.elapsed_nano();  // Record matching time in nanosecond.
        matchTimeList.push_back((double)eventTime / 1000000);
        matchSubList.push_back(matchSubs);

        //cout << matchSubs << endl;
    }
    double avgMatchTime = 0;
    for (int i = 0; i < matchTimeList.size(); i++)avgMatchTime += matchTimeList[i];
    avgMatchTime /= matchTimeList.size();

    cout << "match complete! avgMatchTime = " << avgMatchTime << "ms\n";

    string outputFileName = "match_time.txt";
    ofstream outputfile(outputFileName, ios::out);
    if (!outputfile) {
        cout << "error opening destination file." << endl;
        return 0;
    }
    for (int i = 0; i < pubs; i++) {
        outputfile << matchSubList[i] << " " << matchTimeList[i] << "\n";
    }
    outputfile.close();
    
    a.check(matchSubList);

    return 0;
}