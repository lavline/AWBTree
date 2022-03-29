#include "AWBTree.h"

bool AWBTree::readSubs(string file)
{
    ifstream infile(file, ios::in);
    if (!infile) {
        cout << "error opening source file." << endl;
        return false;
    }
    double tmp;
    for (int i = 0; i < 9; ++i)infile >> tmp;
    subList.resize(subs);
    for (auto& sub : subList) {
        infile >> sub.id;
        infile >> sub.size;
        sub.pred.resize(sub.size);
        for (auto& pred : sub.pred) {
            infile >> pred.att;
            infile >> pred.lowVal;
            infile >> pred.highVal;
        }
    }
    infile.close();
    return true;
}

bool AWBTree::readPus(string file)
{
    ifstream infile(file, ios::in);
    if (!infile) {
        cout << "error opening source file." << endl;
        return false;
    }
    double tmp;
    for (int i = 0; i < 8; ++i)infile >> tmp;
    pubList.resize(pubs);
    for (auto& pub : pubList) {
        infile >> pub.size;
        pub.pairs.resize(pub.size);
        for (auto& pred : pub.pairs) {
            infile >> pred.att;
            infile >> pred.value;
        }
    }
    infile.close();
    return true;
}

void AWBTree::insert(IntervalSub sub)
{
    int id = sub.id;
    counter[id] = sub.size;
    for (auto pred : sub.pred) {
        int pred_att = pred.att;
        int lowVal = pred.lowVal; int highVal = pred.highVal;
        Element ele = { lowVal, highVal, id };
        uint16_t w_id = (highVal - lowVal) / width;
        if (!(w_id > dPoint))++pdr[id];
        sig[pred_att][w_id] = true;
        ltree[pred.att][w_id]->insert(lowVal, ele);
        htree[pred.att][w_id]->insert(highVal, ele);
    }
}

void AWBTree::setbits()
{
    for (int i = 0; i < pdr.size(); i++)
        if (pdr[i] == 0)bit_forward.set(i);
}

void AWBTree::forward(const Pub& pub, int& matchSubs)
{
    vector<uint16_t> c(counter);
    bitset<SUB_NUM> bits;
    for (auto pred : pub.pairs) {
        int pred_att = pred.att; double val = pred.value;
        for (int i = 0; i < width_size; ++i) {
            if (sig[pred_att][i]) {
                double w_min = i * width;
                tree_iter it = ltree[pred_att][i]->upper_bound(val - w_min - width);
                tree_iter sPoint = ltree[pred_att][i]->lower_bound(val - w_min);
                tree_iter end = ltree[pred_att][i]->upper_bound(val);
                for (; it != sPoint; ++it) {
                    for (auto& ele : it.get_value_list())
                        if (!(val > ele.highVal))
                            if (!(--c[ele.subID]))bits.set(ele.subID);
                }
                for (; it != end; ++it) {
                    for (auto& ele : it.get_value_list())
                        if (!(--c[ele.subID]))bits.set(ele.subID);
                }
            }
        }
    }
    matchSubs = bits.count();
}

void AWBTree::forward_o(const Pub& pub, int& matchSubs)
{
    vector<uint16_t> c(counter);
    bitset<SUB_NUM> bits;
    for (auto& pred : pub.pairs) {
        int pred_att = pred.att; double val = pred.value;
        for (int i = 0; i < width_size; ++i) {
            if (sig[pred_att][i]) {
                double w_min = i * width;
                double w_max = w_min + width;
                double v_max = valDom - w_min;
                if (val < valDom) {
                    if (!(val > w_min)) {
                        tree_iter it = ltree[pred_att][i]->begin;
                        tree_iter end = ltree[pred_att][i]->upper_bound(val);
                        for (; it != end; ++it)
                            for (auto& ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits.set(ele.subID);
                    }
                    else if (val < w_max) {
                        tree_iter it = ltree[pred_att][i]->begin;
                        tree_iter sPoint = ltree[pred_att][i]->lower_bound(val - w_min);
                        tree_iter end = ltree[pred_att][i]->upper_bound(val);
                        for (; it != sPoint; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(val > ele.highVal))
                                    if (!(--c[ele.subID]))bits.set(ele.subID);
                        }
                        for (; it != end; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits.set(ele.subID);
                        }
                    }
                    else {
                        tree_iter it = ltree[pred_att][i]->upper_bound(val - w_min - width);
                        tree_iter sPoint = ltree[pred_att][i]->lower_bound(val - w_min);
                        tree_iter end = ltree[pred_att][i]->upper_bound(val);
                        for (; it != sPoint; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(val > ele.highVal))
                                    if (!(--c[ele.subID]))bits.set(ele.subID);
                        }
                        for (; it != end; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits.set(ele.subID);
                        }
                    }
                }
                else {
                    if (!(val < v_max)) {
                        tree_iter it = htree[pred_att][i]->lower_bound(val);
                        tree_iter end = htree[pred_att][i]->end;
                        for (; it != end; ++it)
                            for (auto& ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits.set(ele.subID);
                    }
                    else if (val > (valDom - w_max)) {
                        tree_iter it = ltree[pred_att][i]->lower_bound(val);
                        tree_iter sPoint = ltree[pred_att][i]->upper_bound(val + w_min);
                        tree_iter end = ltree[pred_att][i]->end;
                        for (; it != sPoint; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits.set(ele.subID);
                        }
                        for (; it != end; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(val > ele.highVal))
                                    if (!(--c[ele.subID]))bits.set(ele.subID);
                        }
                    }
                    else {
                        tree_iter it = ltree[pred_att][i]->lower_bound(val);
                        tree_iter sPoint = ltree[pred_att][i]->upper_bound(val + w_min);
                        tree_iter end = ltree[pred_att][i]->lower_bound(val + w_max);
                        for (; it != sPoint; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits.set(ele.subID);
                        }
                        for (; it != end; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(val > ele.highVal))
                                    if (!(--c[ele.subID]))bits.set(ele.subID);
                        }
                    }
                }
            }
        }
    }
    matchSubs = bits.count();
}

void AWBTree::backward(const Pub& pub, int& matchSubs)
{
    bitset<SUB_NUM> bits;
    //Timer start;
    for (auto& pred : pub.pairs) {
        int pred_att = pred.att; double val = pred.value;
        for (int i = 0; i < width_size; ++i) {
            if (sig[pred_att][i]) {
                tree_iter itlow = ltree[pred_att][i]->upper_bound(val);
                tree_iter ithigh = htree[pred_att][i]->lower_bound(val);
                for (; itlow != ltree[pred_att][i]->end; ++itlow)
                    for (auto& ele : itlow.get_value_list())bits.set(ele.subID);
                for (tree_iter it = htree[pred_att][i]->begin; it != ithigh; ++it)
                    for (auto& ele : it.get_value_list())bits.set(ele.subID);
            }
        }
    }
    //int64_t end = start.elapsed_nano();
    //cout << (double)end / 1000000 << endl;
    matchSubs = subs - bits.count();
}

void AWBTree::backward_o(const Pub& pub, int& matchSubs)
{
    bitset<SUB_NUM> bits;
    for (auto& pred : pub.pairs) {
        int pred_att = pred.att; double val = pred.value;
        int i = 0;
        for (; i < width_size >> 1; ++i) {
            if (sig[pred_att][i]) {
                tree_iter itlow = ltree[pred_att][i]->upper_bound(val);
                tree_iter ithigh = htree[pred_att][i]->lower_bound(val);
                for (; itlow != ltree[pred_att][i]->end; ++itlow)
                    for (auto& ele : itlow.get_value_list())bits.set(ele.subID);
                for (tree_iter it = htree[pred_att][i]->begin; it != ithigh; ++it)
                    for (auto& ele : it.get_value_list())bits.set(ele.subID);
            }
        }
        for (; i < width_size; ++i) {
            if (sig[pred_att][i]) {
                double v_min = i * width;
                double v_max = valDom - v_min;
                if (!(val > v_max)) {
                    tree_iter itlow = ltree[pred_att][i]->upper_bound(val);
                    for (; itlow != ltree[pred_att][i]->end; ++itlow)
                        for (auto& ele : itlow.get_value_list())bits.set(ele.subID);
                }
                else if (!(val < v_min)) {
                    tree_iter ithigh = htree[pred_att][i]->lower_bound(val);
                    for (tree_iter it = htree[pred_att][i]->begin; it != ithigh; ++it)
                        for (auto& ele : it.get_value_list())bits.set(ele.subID);
                }
            }
        }
    }
    matchSubs = subs - bits.count();
}

void AWBTree::hybrid(const Pub& pub, int& matchSubs)
{
    vector<uint16_t> c(pdr);
    bitset<SUB_NUM> bits_f(bit_forward);
    bitset<SUB_NUM> bits_b;
    for (auto pred : pub.pairs) {
        int pred_att = pred.att; double val = pred.value;
        int i = 0;
        // forward
        for (; i <= dPoint; ++i) {
            if (sig[pred_att][i]) {
                double w_min = i * width;
                double w_max = w_min + width;
                double v_max = valDom - w_min;
                if (val < valDom) {
                    if (!(val > w_min)) {
                        tree_iter it = ltree[pred_att][i]->begin;
                        tree_iter end = ltree[pred_att][i]->upper_bound(val);
                        for (; it != end; ++it)
                            for (auto ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits_f.set(ele.subID);
                    }
                    else if (val < w_max) {
                        tree_iter it = ltree[pred_att][i]->begin;
                        tree_iter sPoint = ltree[pred_att][i]->lower_bound(val - w_min);
                        tree_iter end = ltree[pred_att][i]->upper_bound(val);
                        for (; it != sPoint; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(val > ele.highVal))
                                    if (!(--c[ele.subID]))bits_f.set(ele.subID);
                        }
                        for (; it != end; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits_f.set(ele.subID);
                        }
                    }
                    else {
                        tree_iter it = ltree[pred_att][i]->upper_bound(val - w_min - width);
                        tree_iter sPoint = ltree[pred_att][i]->lower_bound(val - w_min);
                        tree_iter end = ltree[pred_att][i]->upper_bound(val);
                        for (; it != sPoint; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(val > ele.highVal))
                                    if (!(--c[ele.subID]))bits_f.set(ele.subID);
                        }
                        for (; it != end; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits_f.set(ele.subID);
                        }
                    }
                }
                else {
                    if (!(val < v_max)) {
                        tree_iter it = htree[pred_att][i]->lower_bound(val);
                        tree_iter end = htree[pred_att][i]->end;
                        for (; it != end; ++it)
                            for (auto ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits_f.set(ele.subID);
                    }
                    else if (val > (valDom - w_max)) {
                        tree_iter it = ltree[pred_att][i]->lower_bound(val);
                        tree_iter sPoint = ltree[pred_att][i]->upper_bound(val + w_min);
                        tree_iter end = ltree[pred_att][i]->end;
                        for (; it != sPoint; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits_f.set(ele.subID);
                        }
                        for (; it != end; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(val > ele.highVal))
                                    if (!(--c[ele.subID]))bits_f.set(ele.subID);
                        }
                    }
                    else {
                        tree_iter it = ltree[pred_att][i]->lower_bound(val);
                        tree_iter sPoint = ltree[pred_att][i]->upper_bound(val + w_min);
                        tree_iter end = ltree[pred_att][i]->lower_bound(val + w_max);
                        for (; it != sPoint; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(--c[ele.subID]))bits_f.set(ele.subID);
                        }
                        for (; it != end; ++it) {
                            for (auto& ele : it.get_value_list())
                                if (!(val > ele.highVal))
                                    if (!(--c[ele.subID]))bits_f.set(ele.subID);
                        }
                    }
                }
            }
        }
        // backward
        for (; i < width_size >> 1; ++i) {
            if (sig[pred_att][i]) {
                tree_iter itlow = ltree[pred_att][i]->upper_bound(val);
                tree_iter ithigh = htree[pred_att][i]->lower_bound(val);
                for (; itlow != ltree[pred_att][i]->end; ++itlow)
                    for (auto ele : itlow.get_value_list())bits_b.set(ele.subID);
                for (tree_iter it = htree[pred_att][i]->begin; it != ithigh; ++it)
                    for (auto ele : it.get_value_list())bits_b.set(ele.subID);
            }
        }
        for (; i < width_size; ++i) {
            if (sig[pred_att][i]) {
                double v_min = i * width;
                double v_max = valDom - v_min;
                if (!(val > v_max)) {
                    tree_iter itlow = ltree[pred_att][i]->upper_bound(val);
                    for (; itlow != ltree[pred_att][i]->end; ++itlow)
                        for (auto ele : itlow.get_value_list())bits_b.set(ele.subID);
                }
                else if (!(val < v_min)) {
                    tree_iter ithigh = htree[pred_att][i]->lower_bound(val);
                    for (tree_iter it = htree[pred_att][i]->begin; it != ithigh; ++it)
                        for (auto ele : it.get_value_list())bits_b.set(ele.subID);
                }
            }
        }
    }
    matchSubs = (bits_f & bits_b.flip()).count();
}

void AWBTree::hybrid_a(const Pub& pub, int& matchSubs)
{
}

size_t AWBTree::memory()
{
    size_t mem = sizeof(bptree*) * dim * width_size * 2;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < width_size; j++) {
            mem += ltree[i][j]->memory();
            mem += htree[i][j]->memory();
        }
    }
    return mem;
}

bool AWBTree::check(vector<uint32_t>& matchSubList)
{
    cout << "check start...\n";
    int i = 0;
    for (auto& pub : pubList) {
        int n = 0;
        vector<int> index(dim, 0);
        for (int j = 0; j < pub.pairs.size(); j++)index[pub.pairs[j].att] = j + 1;
        for (auto& sub : subList) {
            int match = 0;
            for (auto& spred : sub.pred) {
                int pid = index[spred.att];
                if (index[spred.att] == 0)continue;
                if (pub.pairs[--pid].value < spred.lowVal || pub.pairs[pid].value > spred.highVal)break;
                ++match;
            }
            if (match == sub.size)++n;
        }
        if (n != matchSubList[i]) {
            cout << "pub " << i << " match error. " << "True match is " << n << ", but result is " << matchSubList[i] << endl;
            return false;
        }
        ++i;
    }
    cout << "result is correct.\n";
    return true;
}
