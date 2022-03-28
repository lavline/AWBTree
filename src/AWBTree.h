#pragma once
#include <bitset>
#include <fstream>
#include "BplusTree.h"
#include "data_structure.h"
#include "ThreadPool.h"
#define SUB_NUM 1000000

using bptree = BplusTree::bplus_tree<double, Element>;
using bptree_node = BplusTree::Node<double, Element>;
using tree_iter = BplusTree::bplus_tree<double, Element>::iterator;

class AWBTree {
private:
	vector<vector<bptree*>> ltree;
	vector<vector<bptree*>> htree;

	uint32_t dim;  // 维度
	uint32_t atts; // 属性数
	uint16_t width_size; // 宽度单元数目
	uint16_t branch; // btree最大分支数
	vector<uint16_t> pdr;  // 记录正向匹配需要匹配的谓词数目
	vector<vector<bool>> sig;
	vector<uint16_t> counter;

	int valDom; // 值域空间
	double width;  // 宽度单元大小
	uint16_t dPoint; // 任务切分点
	uint32_t subs; //订阅数
	uint32_t pubs; //

public:
	vector<IntervalSub> subList;
	vector<Pub> pubList;

public:
	AWBTree(uint32_t _dim, uint32_t _atts, uint32_t _subs, uint32_t _pubs, uint16_t _width_size, uint16_t _branch, int _valDom, uint16_t _dPoint)
		:dim(_dim), atts(_atts), subs(_subs), pubs(_pubs), width_size(_width_size), branch(_branch), width(_valDom / _width_size),
		valDom(_valDom), sig(_dim, vector<bool>(_width_size, 0)), counter(_subs, 0), dPoint(_dPoint* _width_size)
	{
		ltree.resize(dim); htree.resize(dim);
		for (int i = 0; i < dim; i++) {
			ltree[i].resize(width_size); htree[i].resize(width_size);
			for (int j = 0; j < width_size; j++) {
				ltree[i][j] = new bptree(branch, new bptree_node);
				htree[i][j] = new bptree(branch, new bptree_node);
			}
		}
	}

	void readSubs(string file);
	void readPus(string file);

	void insert(IntervalSub sub);
	void erase(IntervalSub sub);

	void forward(const Pub& pub, int& matchSubs);
	void forward_o(const Pub& pub, int& matchSubs);
	void backward(const Pub& pub, int& matchSubs);
	void backward_o(const Pub& pub, int& matchSubs);
	void hybrid(const Pub& pub, int& matchSubs);
	void hybrid_a(const Pub& pub, int& matchSubs);

	size_t memory();
	bool check(vector<uint32_t>& matchSubList);
};