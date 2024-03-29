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
        //Element ele = { lowVal, highVal, id };
        lowTreeEle ele = { highVal, id };
        uint16_t w_id = (highVal - lowVal) / width;
        if (!(w_id > dPoint))++pdr[id];
        sig[pred_att][w_id] = true;
        ltree[pred.att][w_id]->insert(lowVal, ele);
        htree[pred.att][w_id]->insert(highVal, id);
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
    for (auto& pred : pub.pairs) {
        int pred_att = pred.att; double val = pred.value;
        for (int i = 0; i < width_size; ++i) {
            if (sig[pred_att][i]) {
                double w_min = i * width;
                ltree_iter it = ltree[pred_att][i]->upper_bound(val - w_min - width);
                ltree_iter sPoint = ltree[pred_att][i]->lower_bound(val - w_min);
                ltree_iter end = ltree[pred_att][i]->upper_bound(val);
                for (; it != sPoint; ++it) {
                    lowTreeEle ele = it.get_val();
                    if (!(val > ele.highVal))
                        if (!(--c[ele.subID]))bits[ele.subID] = true;
                }
                for (; it != end; ++it) {
                    lowTreeEle ele = it.get_val();
                    if (!(--c[ele.subID]))bits[ele.subID] = true;
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
    //Timer start;
    for (auto& pred : pub.pairs) {
        int pred_att = pred.att; int val = pred.value;
        int i = 0;
        for (; i < width_size >> 1; ++i) {
            if (sig[pred_att][i]) {
                double w_min = i * width;
                ltree_iter it = ltree[pred_att][i]->upper_bound(val - w_min - width);
                ltree_iter sPoint = ltree[pred_att][i]->lower_bound(val - w_min);
                ltree_iter end = val < (valDom - w_min) ? ltree[pred_att][i]->upper_bound(val) : ltree[pred_att][i]->end();
                for (; it != sPoint; ++it) {
                    //Element ele = it.get_value();
                    lowTreeEle ele = it.get_val();
                    if (!(val > ele.highVal))
                        if (!(--c[ele.subID]))bits[ele.subID] = true;
                }
                for (; it != end; ++it) {
                    //Element ele = it.get_value();
                    lowTreeEle ele = it.get_val();
                    if (!(--c[ele.subID]))bits[ele.subID] = true;
                }
            }
        }
        for (; i < width_size; ++i) {
            if (sig[pred_att][i]) {
                double w_min = i * width;
                double w_max = w_min + width;
                double v_max = valDom - w_min;
                if (val > w_min) {
                    ltree_iter it = ltree[pred_att][i]->upper_bound(val - w_min - width);
                    ltree_iter sPoint = ltree[pred_att][i]->lower_bound(val - w_min);
                    ltree_iter end = ltree[pred_att][i]->end();
                    for (; it != sPoint; ++it) {
                        lowTreeEle ele = it.get_val();
                        if (!(val > ele.highVal))
                            if (!(--c[ele.subID]))bits[ele.subID] = true;
                    }
                    for (; it != end; ++it) {
                        lowTreeEle ele = it.get_val();
                        if (!(--c[ele.subID]))bits[ele.subID] = true;
                    }
                }
                else if (val > v_max) {
                    ltree_iter it = ltree[pred_att][i]->begin();
                    ltree_iter end = ltree[pred_att][i]->end();
                    for (; it != end; ++it) {
                        lowTreeEle ele = it.get_val();
                        if (!(--c[ele.subID]))bits[ele.subID] = true;
                    }
                }
                else {
                    ltree_iter it = ltree[pred_att][i]->begin();
                    ltree_iter end = ltree[pred_att][i]->upper_bound(val);
                    for (; it != end; ++it) {
                        lowTreeEle ele = it.get_val();
                        if (!(--c[ele.subID]))bits[ele.subID] = true;
                    }
                }
                //if (!(val > w_min)) {
                //    ltree_iter it = ltree[pred_att][i]->begin();
                //    ltree_iter end = ltree[pred_att][i]->upper_bound(val);
                //    for (; it != end; ++it) {
                //        //Element ele = it.get_value();
                //        lowTreeEle ele = it.get_val();
                //        if (!(--c[ele.subID]))bits[ele.subID] = true;
                //    }
                //}
                //else if (val < w_max) {
                //    ltree_iter it = ltree[pred_att][i]->begin();
                //    ltree_iter sPoint = ltree[pred_att][i]->lower_bound(val - w_min);
                //    ltree_iter end = ltree[pred_att][i]->upper_bound(val);
                //    for (; it != sPoint; ++it) {
                //        //Element ele = it.get_value();
                //        lowTreeEle ele = it.get_val();
                //        if (!(val > ele.highVal))
                //            if (!(--c[ele.subID]))bits[ele.subID] = true;
                //    }
                //    for (; it != end; ++it) {
                //        //Element ele = it.get_value();
                //        lowTreeEle ele = it.get_val();
                //        if (!(--c[ele.subID]))bits[ele.subID] = true;
                //    }
                //}
                //else {
                //    ltree_iter it = ltree[pred_att][i]->upper_bound(val - w_min - width);
                //    ltree_iter sPoint = ltree[pred_att][i]->lower_bound(val - w_min);
                //    ltree_iter end = ltree[pred_att][i]->upper_bound(val);
                //    for (; it != sPoint; ++it) {
                //        //Element ele = it.get_value();
                //        lowTreeEle ele = it.get_val();
                //        if (!(val > ele.highVal))
                //            if (!(--c[ele.subID]))bits[ele.subID] = true;
                //    }
                //    for (; it != end; ++it) {
                //        //Element ele = it.get_value();
                //        lowTreeEle ele = it.get_val();
                //        if (!(--c[ele.subID]))bits[ele.subID] = true;
                //    }
                //}
            }
        }
    }
    //int64_t end = start.elapsed_nano();
    //cout << (double)end / 1000000 << " ms\n";
    matchSubs = bits.count();
}

void AWBTree::backward(const Pub& pub, int& matchSubs)
{
    bitset<SUB_NUM> bits;
    //Timer start;
    for (auto& pred : pub.pairs) {
        int pred_att = pred.att; int val = pred.value;
        for (int i = 0; i < width_size; ++i) {
            if (sig[pred_att][i]) {
                ltree_iter itlow = ltree[pred_att][i]->upper_bound(val);
                htree_iter ithigh = htree[pred_att][i]->lower_bound(val);
                for (; itlow != ltree[pred_att][i]->end(); ++itlow) {
                    //bits.set(itlow.get_val().subID);
                    bits[itlow.get_val().subID] = true;
                }
                for (htree_iter it = htree[pred_att][i]->begin(); it != ithigh; ++it) {
                    //bits.set(it.get_val());
                    bits[it.get_val()] = true;
                }
            }
        }
    }
    //int64_t end = start.elapsed_nano();
    //cout << (double)end / 1000000 << " ms\n";
    //Timer start;
    matchSubs = subs - bits.count();
    //int64_t end = start.elapsed_nano();
    //cout << (double)end / 1000 << "um\n";
}

void AWBTree::backward_o(const Pub& pub, int& matchSubs)
{
    bitset<SUB_NUM> bits;
    for (auto& pred : pub.pairs) {
        int pred_att = pred.att; int val = pred.value;
        int i = 0;
        for (; i < width_size >> 1; ++i) {
            if (sig[pred_att][i]) {
                ltree_iter itlow = ltree[pred_att][i]->upper_bound(val);
                htree_iter ithigh = htree[pred_att][i]->lower_bound(val);
                for (; itlow != ltree[pred_att][i]->end(); ++itlow)
                    bits[itlow.get_val().subID] = true;
                for (htree_iter it = htree[pred_att][i]->begin(); it != ithigh; ++it)
                    bits[it.get_val()] = true;
            }
        }
        for (; i < width_size; ++i) {
            if (sig[pred_att][i]) {
                double v_min = i * width;
                double v_max = valDom - v_min;
                if (!(val > v_max)) {
                    ltree_iter itlow = ltree[pred_att][i]->upper_bound(val);
                    for (; itlow != ltree[pred_att][i]->end(); ++itlow)
                        bits[itlow.get_val().subID] = true;
                }
                else if (!(val < v_min)) {
                    htree_iter ithigh = htree[pred_att][i]->lower_bound(val);
                    for (htree_iter it = htree[pred_att][i]->begin(); it != ithigh; ++it)
                        bits[it.get_val()] = true;
                }
            }
        }
    }
    matchSubs = subs - bits.count();
}

void AWBTree::hybrid(const Pub& pub, int& matchSubs)
{
    //Timer start;
    vector<uint16_t> c(pdr);
    bitset<SUB_NUM> bits_f(bit_forward);
    bitset<SUB_NUM> bits_b;
    //int64_t end = start.elapsed_nano();
    //cout << (double)end / 1000 << "um\n";
    for (auto& pred : pub.pairs) {
        int pred_att = pred.att; int val = pred.value;
        int i = 0;
        // forward
        for (; i <= dPoint; ++i) {
            if (sig[pred_att][i]) {
                double w_min = i * width;
                ltree_iter it = ltree[pred_att][i]->upper_bound(val - w_min - width);
                ltree_iter sPoint = ltree[pred_att][i]->lower_bound(val - w_min);
                ltree_iter end = val < (valDom - w_min) ? ltree[pred_att][i]->upper_bound(val) : ltree[pred_att][i]->end();
                for (; it != sPoint; ++it) {
                    lowTreeEle ele = it.get_val();
                    if (!(val > ele.highVal))
                        if (!(--c[ele.subID]))bits_f[ele.subID] = true;
                }
                for (; it != end; ++it) {
                    lowTreeEle ele = it.get_val();
                    if (!(--c[ele.subID]))bits_f[ele.subID] = true;
                }
            }
        }
        // backward
        for (; i < width_size >> 1; ++i) {
            if (sig[pred_att][i]) {
                ltree_iter itlow = ltree[pred_att][i]->upper_bound(val);
                htree_iter ithigh = htree[pred_att][i]->lower_bound(val);
                for (; itlow != ltree[pred_att][i]->end(); ++itlow)
                    bits_b[itlow.get_val().subID] = true;
                for (htree_iter it = htree[pred_att][i]->begin(); it != ithigh; ++it)
                    bits_b[it.get_val()] = true;
            }
        }
        for (; i < width_size; ++i) {
            if (sig[pred_att][i]) {
                double v_min = i * width;
                double v_max = valDom - v_min;
                if (!(val > v_max)) {
                    ltree_iter itlow = ltree[pred_att][i]->upper_bound(val);
                    for (; itlow != ltree[pred_att][i]->end(); ++itlow)
                        bits_b[itlow.get_val().subID] = true;
                }
                else if (!(val < v_min)) {
                    htree_iter ithigh = htree[pred_att][i]->lower_bound(val);
                    for (htree_iter it = htree[pred_att][i]->begin(); it != ithigh; ++it)
                        bits_b[it.get_val()] = true;
                }
            }
        }
    }
    //Timer start;
    matchSubs = (bits_f & bits_b.flip()).count();
    /*int64_t end = start.elapsed_nano();
    cout << (double)end / 1000 << "um\n";*/
}

void AWBTree::hybrid_a(const Pub& pub, int& matchSubs)
{
}

size_t AWBTree::memory()
{
    size_t mem = sizeof(lbptree*) * dim * width_size * 2;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < width_size; j++) {
            mem += ltree[i][j]->mem_size();
            mem += htree[i][j]->mem_size();
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
