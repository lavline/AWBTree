#pragma once
//#ifndef BPLUSTREE_H
//#define BPLUSTREE_H
#include <vector>
#include <stdint.h>
#include <iostream>

namespace BplusTree {
	template <class key_type, class value_type>
	class Node {
	private:
		template <class key, class value>
		friend class bplus_tree;
		bool NON_LEAF;
		std::vector<key_type> key;  // key
		std::vector<std::vector<value_type>> value;  // value: each value type have a value list to store the value with same key
		std::vector<Node<key_type, value_type>*> child;
		//class Node<key_type, value_type>* parent;
		Node<key_type, value_type>* next;
	public:
		Node() : next(nullptr), NON_LEAF(false) {
			//parent = nullptr;
		}
		Node(const std::vector<key_type>& _key, const std::vector<std::vector<value_type>>& _value, const int mid)
			:key(_key.begin() + mid + 1, _key.end()), value(_value.begin() + mid + 1, _value.end()), next(nullptr), NON_LEAF(false) {}
		Node(const std::vector<key_type>& _key, const std::vector<Node<key_type, value_type>*>& _child, const int mid)
			:key(_key.begin() + mid + 1, _key.end()), child(_child.begin() + mid + 1, _child.end()), next(nullptr), NON_LEAF(true) {}
	};

	template <class key_type, class value_type>
	class bplus_tree {
	public:
		class iterator;
		const iterator begin;
		const iterator end;

	private:
		const int branch;
		const int min_branch;
		Node<key_type, value_type>* root;

	public:
		bplus_tree(const uint16_t _branch, Node<key_type, value_type>* _root) 
			:branch(_branch), min_branch((branch - 1) >> 1), root(_root), begin(_root, 0), end() {}
		~bplus_tree() {
			clear(root);
		}
	private:
		size_t nodes_memory(Node<key_type, value_type>* node) {
			size_t mem = sizeof(Node<key_type, value_type>) + node->key.size() * sizeof(key_type)
				+ node->child.size() * sizeof(Node<key_type, value_type>*) + node->value.size() * sizeof(std::vector<value_type>);
			for (int i = 0; i < node->value.size(); i++)mem += node->value[i].size() * sizeof(value_type);
			for (int i = 0; i < node->child.size(); i++)mem += nodes_memory(node->child[i]);
			return mem;
		}
	public:
		size_t memory() {
			size_t mem = sizeof(bplus_tree<key_type, value_type>);
			Node<key_type, value_type>* node = root;
			return mem + nodes_memory(root);
		}
		void insert(const key_type key, const value_type& value) {
			// insert key to leaf node
			Node<key_type, value_type>* node = root;
			std::vector<Node<key_type, value_type>*> traverse_node;
			std::vector<int> traverse_index;
			int i = 0;
			size_t node_size = node->key.size();
			// find leaf node
			while (node->NON_LEAF) {
				while (i < node_size && key > node->key[i])++i;
				traverse_node.emplace_back(node);
				traverse_index.emplace_back(i);
				node = node->child[i];
				node_size = node->key.size();
				i = 0;
			}
			// find key
			while (i < node_size && key > node->key[i])++i;
			// if not exist the key, build a new <key, value> pair
			if (node_size == 0|| i == node_size || key != node->key[i]) {
				node->key.insert(node->key.begin() + i, key);
				node->value.insert(node->value.begin() + i, std::vector<value_type>(1, value));
				// node split from leaf node if need
				int mid = min_branch;
				if (node_size + 1 == branch) {
					key_type mid_key = node->key[mid];
					// create right leaf node
					Node<key_type, value_type>* new_node = new Node<key_type, value_type>(node->key, node->value, mid);
					//new_node->value.insert(new_node->value.end(), node->value.begin() + mid + 1, node->value.end());
					new_node->next = node->next;
					node->next = new_node;
					// resize left node
					node->key.resize(mid + 1);
					node->value.resize(mid + 1);
					// create new root
					if (traverse_node.size() == 0) {
						Node<key_type, value_type>* new_root = new Node<key_type, value_type>;
						new_root->key.emplace_back(mid_key);
						new_root->child.emplace_back(node);
						new_root->child.emplace_back(new_node);
						new_root->NON_LEAF = true;
						root = new_root;
						return;
					}
					// else modify parent node
					Node<key_type, value_type>* parent = traverse_node.back();
					traverse_node.pop_back();
					int idx = traverse_index.back();
					traverse_index.pop_back();
					parent->key.insert(parent->key.begin() + idx, mid_key);
					parent->child.insert(parent->child.begin() + idx + 1, new_node);
					node = parent;
					node_size = node->key.size();
					// internal node split if need
					while (node_size == branch) {
						key_type mid_key = node->key[mid];
						// creat right internal node
						Node<key_type, value_type>* new_node = new Node<key_type, value_type>(node->key, node->child, mid);
						// resize left interal node
						node->key.resize(mid);
						node->child.resize(mid + 1);
						// creat new root
						if (traverse_node.size() == 0) {
							Node<key_type, value_type>* new_root = new Node<key_type, value_type>;
							new_root->key.emplace_back(mid_key);
							new_root->child.emplace_back(node);
							new_root->child.emplace_back(new_node);
							new_root->NON_LEAF = true;
							root = new_root;
							return;
						}
						// else modify parent node
						Node<key_type, value_type>* parent = traverse_node.back();
						traverse_node.pop_back();
						int idx = traverse_index.back();
						traverse_index.pop_back();
						parent->key.insert(parent->key.begin() + idx, mid_key);
						parent->child.insert(parent->child.begin() + idx + 1, new_node);
						node = parent;
						node_size = node->key.size();
					}
				}
			}
			else node->value[i].emplace_back(value);
		}
		void erase(){}

		iterator lower_bound(const key_type key) const {
			iterator it;
			Node<key_type, value_type>* node = root;
			size_t node_size = node->key.size();
			// find leaf node
			while (node->NON_LEAF) {
				int i = 0;
				while (i < node_size && key > node->key[i])++i;
				node = node->child[i];
				node_size = node->key.size();
			}
			// find the first key >= the given key
			for (int i = 0; i < node_size; ++i) {
				if (key > node->key[i])continue;
				it.node = node;
				it.index = i;
				return it;
			}
			return this->end;
		}
		iterator upper_bound(const key_type key) const {
			iterator it = lower_bound(key);
			if (it != this->end) ++it;
			return it;
		}

		void clear(Node<key_type,value_type>* node){
			for (int i = 0; i < node->child.size(); ++i) {
				clear(node->child[i]);
			}
			delete node;
		}

		class iterator {
		private:
			friend class bplus_tree;
			Node<key_type, value_type>* node;
			int index;

		public:
			iterator():node(nullptr),index(0) {}
			iterator(Node<key_type, value_type>* _node, int _index) :node(_node), index(_index){}

			key_type& get_key() const {
				if (node == nullptr)
					throw std::out_of_range("B+Tree: iterator is out of range");
				return this->node->key[index];
			}
			std::vector<value_type>& get_value_list() const {
				if (node == nullptr)
					throw std::out_of_range("B+Tree: iterator is out of range");
				return this->node->value[index];
			}

			// prefix increment operator (++it). it returns a referece to the incremented iterator, avoiding the copy.
			const iterator& operator++(){
				if (node == nullptr)
					throw std::out_of_range("B+Tree: iterator is out of range");
				if (index + 1 < node->key.size())
					index++;
				else{
					index = 0;
					node = node->next;
				}
				return *this;
			}
			bool operator==(const iterator& it) const {
				return (this->node == it.node && this->index == it.index);
			}
			bool operator!=(const iterator& it) const {
				return (this->node != it.node || this->index != it.index);
			}
		};
	};
}

//#endif // BPLUSTREE_H
