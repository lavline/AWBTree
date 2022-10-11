#pragma once
#include <bitset>
#include <fstream>
#include <mutex>
#include <atomic>
#include "BplusTree.h"
#include "chrono_time.h"
#include "data_structure.h"
#include "ThreadPool.h"
#define SUB_NUM 1000000


//using lbptree = BplusTree::bplus_tree<int, lowTreeEle>;
//using lbptree_node = BplusTree::Node<int, lowTreeEle>;
//using ltree_iter = BplusTree::bplus_tree<int, lowTreeEle>::iterator;
//using hbptree = BplusTree::bplus_tree<int, int>;
//using hbptree_node = BplusTree::Node<int, int>;
//using htree_iter = BplusTree::bplus_tree<int, int>::iterator;
using lbptree = BplusTree::Tree<int, lowTreeEle>;
using lbptree_node = BplusTree::Node<int, lowTreeEle>;
using ltree_iter = BplusTree::Tree<int, lowTreeEle>::iterator;
using hbptree = BplusTree::Tree<int, int>;
using hbptree_node = BplusTree::Node<int, int>;
using htree_iter = BplusTree::Tree<int, int>::iterator;

class AWBTree {
private:
	vector<vector<lbptree*>> ltree;
	vector<vector<hbptree*>> htree;

	uint32_t dim;  // ά��
	uint32_t atts; // ������
	uint16_t width_size; // ��ȵ�Ԫ��Ŀ
	uint16_t branch; // btree����֧��
	vector<uint16_t> pdr;  // ��¼����ƥ����Ҫƥ���ν����Ŀ
	vector<vector<bool>> sig;
	vector<uint16_t> counter;
	//vector<atomic_uint16_t> atomic_counter;
	bitset<SUB_NUM> bit_forward;

	int valDom; // ֵ��ռ�
	double width;  // ��ȵ�Ԫ��С
	uint16_t dPoint; // �����зֵ�
	uint32_t subs; //������
	uint32_t pubs; //

public:
	vector<IntervalSub> subList;
	vector<Pub> pubList;

public:
	AWBTree(uint32_t _dim, uint32_t _atts, uint32_t _subs, uint32_t _pubs, uint16_t _width_size, uint16_t _branch, int _valDom, float _dPoint)
		:dim(_dim), atts(_atts), subs(_subs), pubs(_pubs), width_size(_width_size), branch(_branch), width(_valDom / _width_size),
		valDom(_valDom), sig(_dim, vector<bool>(_width_size, 0)), counter(_subs, 0), dPoint(_dPoint * _width_size), pdr(_subs, 0)
	{
		ltree.resize(dim); htree.resize(dim);
		for (int i = 0; i < dim; i++) {
			ltree[i].resize(width_size); htree[i].resize(width_size);
			for (int j = 0; j < width_size; j++) {
				//ltree[i][j] = new lbptree(branch, new lbptree_node);
				//htree[i][j] = new hbptree(branch, new hbptree_node);
				ltree[i][j] = new lbptree(branch);
				htree[i][j] = new hbptree(branch);
			}
		}
	}

	bool readSubs(string file);
	bool readPus(string file);

	void insert(IntervalSub sub);
	void setbits();
	void erase(IntervalSub sub);

	void forward(const Pub& pub, int& matchSubs);
	void forward_o(const Pub& pub, int& matchSubs);
	void forward_p(const Pub& pub, int& matchSubs, ThreadPool& pool, int pdegree);
	void backward(const Pub& pub, int& matchSubs);
	void backward_o(const Pub& pub, int& matchSubs);
	void backward_p(const Pub& pub, int& matchSubs, ThreadPool& pool, int pdegree);
	void hybrid(const Pub& pub, int& matchSubs);
	void hybrid_a(const Pub& pub, int& matchSubs);
	void hybrid_p(const Pub& pub, int& matchSubs, ThreadPool& pool, int pdegree);

	size_t memory();
	bool check(vector<uint32_t>& matchSubList);
};