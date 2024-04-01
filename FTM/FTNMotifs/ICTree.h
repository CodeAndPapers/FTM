#pragma once
#include "stdafx.h"
#include "Interval.h"
#include "Util.h"

/*IC trees*/
using Intv = Interval<int, Label>;

class ICTNode {
public:

	ICTNode(const ICTNode& node)
		:key(node.key), left(node.left), right(node.right), parent(node.parent),
		color(node.color) {}

	ICTNode(const Intv& val)
		:key(val), left(NULL), right(NULL), parent(NULL), color(true) {}

	ICTNode(const Intv& val, ICTNode* &parent)
		:key(val), left(NULL), right(NULL), parent(parent), color(true) {}

	friend ostream& operator<<(ostream& out, const ICTNode*node) {
		out << "key:" << node->key
			<< "\tcolor:" << node->color;
		return out;
	}

	~ICTNode() = default;

	ICTNode*left, *right, *parent;//left\right child, parent node  24byte
	/*range : max right endpoint of interval in the subtree and
	min left endpoint of interval in the subtree*/
	//int maxV, minV; // 8byte
	Intv key;//interval // 12byte
	bool color;//true(RED) or false(BLACK) 1byte
};

class ICTree {
public:
	ICTree() :root(NULL) {}

	ICTree(const ICTree& tree) {
		root = DBG_NEW ICTNode(tree.root->key);
		root->color = tree.root->color;
		if (tree.root->left) {
			copy(root, root->left, tree.root->left);
		}
		if (tree.root->right) {
			copy(root, root->right, tree.root->right);
		}
	}

	/*search if one node contain the time
	and return a interval if true 
	(used in updating IC trees)*/
	bool search(Intv& intv, int time);

	/*update interval's endT of nodes
	(used in updating IC trees)*/
	void updateNode(Intv& intv, int newEndT);

	/*encapsulation of insertion 
	(used in creating IC trees)*/
	void insert(const Intv& val);

	/*search if there exists a interval containing [startT,endT],
	and return the maximum interval [intvStartT, intvEndT] containing 
	[startT,endT] where the edge keep the same label edgeType*/
	bool containQuery(const int startT, const int endT, 
		Label& edgeType, int& intvStartT, int& intvEndT);

	void print();

	~ICTree() {
		release(root);
	}

	/*copy the information in copynode to node*/
	void copy(ICTNode*& parent,
		ICTNode*& node, ICTNode*& copynode);


private:
	ICTNode* root;

	#pragma region maintain interval tree
	//set root node
	inline void setRoot(Intv& val) {
		this->root = DBG_NEW ICTNode(val);
	}

	//rotate right
	void rRotate(ICTNode* &node);

	//rotate left
	void lRotate(ICTNode* &node);

	//maintain the property of red-black
	void insertFixUp(ICTNode* &node);

	//release memory
	void release(ICTNode* &node);
	#pragma endregion
};