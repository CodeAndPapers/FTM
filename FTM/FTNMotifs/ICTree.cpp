#include "ICTree.h"
#include "stdafx.h"

#pragma region maintain interval tree
//rotate right
void ICTree::rRotate(ICTNode* &mynode) {
	ICTNode* node = mynode;
	if (!node)return;
	ICTNode* temp = node->left;
	if (!temp)return;
	node->left = temp->right;//right subtree of temp put in root->left
	if (temp->right)temp->right->parent = node;//the right node!=NULL
	temp->right = node;//root put in right subtree of temp
	if (node->parent&&node->parent->left == node) {
		node->parent->left = temp;
	}
	else if (node->parent) {
		node->parent->right = temp;
	}
	temp->parent = node->parent;
	node->parent = temp;
}

//rotate left
void ICTree::lRotate(ICTNode* &mynode) {
	ICTNode* node = mynode;
	if (!node)return;
	ICTNode* temp = node->right;
	if (!temp)return;
	node->right = temp->left;//left subtree of temp put in root->right
	if (temp->left)temp->left->parent = node;//the right node!=NULL
	temp->left = node;//root put in left subtree of temp
	if (node->parent&&node->parent->left == node) {
		node->parent->left = temp;
	}
	else if (node->parent) {
		node->parent->right = temp;
	}
	temp->parent = node->parent;
	node->parent = temp;
}

//maintain the property of red-black
void ICTree::insertFixUp(ICTNode* &mynode) {
	ICTNode* node = mynode;
	while (node->parent != NULL && node->parent->color) {
		if (node->parent == node->parent->parent->left) {//left subtree
			/*if node->parent exists, node->parent->parent must exist,
			because node->parent->color==red and
			node->parent can't be the root
			*/
			ICTNode* uncle = node->parent->parent->right;
			if (uncle&&uncle->color) {//uncle can be NULL
				uncle->color = false;
				node->parent->color = false;
				node->parent->parent->color = true;
				node = node->parent->parent;
			}
			else {
				if (node == node->parent->right) {//right subtree
					node = node->parent;
					lRotate(node);
				}
				node->parent->color = false;
				node->parent->parent->color = true;
				rRotate(node->parent->parent);
			}
		}
		else {//right subtree
			ICTNode* uncle = node->parent->parent->left;
			if (uncle&&uncle->color == true) {//uncle can be NULL
				uncle->color = false;
				node->parent->color = false;
				node->parent->parent->color = true;
				node = node->parent->parent;
			}
			else {
				if (node == node->parent->left) {//right subtree
					node = node->parent;
					rRotate(node);
				}
				node->parent->color = false;
				node->parent->parent->color = true;
				lRotate(node->parent->parent);
			}
		}
	}
	if (this->root->parent)//change the root ptr
		this->root = this->root->parent;
	this->root->color = false;//set the color of root
}

void ICTree::release(ICTNode* &node) {
	if (node) {
		release(node->right);
		release(node->left);
		delete node;
		node = NULL;
	}
}
#pragma endregion

/*encapsulation of insertion (used in creating IC trees)*/
void ICTree::insert(const Intv& val) {
	if (!this->root) {
		this->root = DBG_NEW ICTNode(val);
		insertFixUp(this->root);
		return;
	}
	ICTNode *now = this->root;
	ICTNode *node;
	while (now != NULL) {
		if (val > now->key) {//right subtree
			if (now->right) {
				now = now->right;
			}
			else {//null
				now->right = DBG_NEW ICTNode(val, now);
				node = now->right;
				break;
			}
		}
		else {//left subtree
			if (now->left) {
				now = now->left;
			}
			else {//null
				now->left = DBG_NEW ICTNode(val, now);
				node = now->left;
				break;
			}
		}
	}
	insertFixUp(node);
}

/*search if there exists a interval containing [startT,endT],
and return the maximum interval [intvStartT, intvEndT] containing
[startT,endT] where the edge keep the same label edgeType*/
bool ICTree::containQuery(const int startT, const int endT,
	Label& edgeType, int& intvStartT, int& intvEndT){
	ICTNode *now = this->root;
	if (!now)return false;
	queue<ICTNode*> q;
	q.push(now);
	Intv* intvPtr;
	int tempStartT, tempEndT;
	while (!q.empty()) {
		now = q.front();
		intvPtr = &now->key;
		tempStartT = intvPtr->startT;
		tempEndT = intvPtr->endT;
		if (startT <= tempEndT &&tempStartT <= endT) {//intersection relationship
			if (tempStartT <= startT && endT <= tempEndT) {//contain relationship
				edgeType = intvPtr->value;
				intvStartT = tempStartT;
				intvEndT = tempEndT;
				return true;
			}
			edgeType = intvPtr->value;
			intvStartT = tempStartT;
			intvEndT = tempEndT;
			return false;
		}
		q.pop();
		if (now->left&&endT <= intvPtr->startT)//left subtree
			q.push(now->left);
		else if (now->right)//right subtree
			q.push(now->right);
	}
	return false;
}

void ICTree::print() {
	if (!this->root) { cout << "empty tree" << endl; return; }
	queue<ICTNode* > q;
	q.push(this->root);
	while (!q.empty()) {
		ICTNode* now = q.front();
		q.pop();
		cout << now << endl;
		if (now->left)q.push(now->left);
		if (now->right)q.push(now->right);
	}
}

/*copy the information in copynode to node*/
void ICTree::copy(ICTNode* &parent,
	ICTNode* &node, ICTNode* &copynode) {
	node->color = copynode->color;
	node->parent = parent;
	if (!copynode->left) {
		node->left = NULL;
	}
	else {
		node->left = DBG_NEW ICTNode(copynode->left->key);
		copy(node, node->left, copynode->left);
	}
	if (!copynode->right) {
		node->right = NULL;
	}
	else {
		node->right = DBG_NEW ICTNode(copynode->right->key);
		copy(node, node->right, copynode->right);
	}
}

/*search if one node contain the time
and return a interval if true
(used in updating IC trees)*/
bool ICTree::search(Intv& intv, int time) {
	ICTNode *now = this->root;
	if (!now) return false;
	queue<ICTNode*> q;
	q.push(now);
	Intv* intvPtr;
	while (!q.empty()) {
		now = q.front();
		q.pop();
		intvPtr = &now->key;
		if (time <= intvPtr->endT&&intvPtr->startT <= time) {
			intv.startT = intvPtr->startT;
			intv.endT = intvPtr->endT;
			intv.value = intvPtr->value;
			return true;
		}
		if (now->left&&time < intvPtr->startT)//left subtree
			q.push(now->left);
		else if (now->right)//right subtree
			q.push(now->right);
	}
	return false;
}

/*update interval's endT of nodes
	(used in updating IC trees)*/
void ICTree::updateNode(Intv& intv, int newEndT) {
	ICTNode *now = this->root;
	while (now) {
		if (intv.startT > now->key.startT) {//right subtree
			if (now->right) {
				now = now->right;
			}
			else {//impossible condition
				cout << "interval doesn't exist" << endl;
				exit(-1);
			}
		}
		else if (intv.startT < now->key.startT) {//left subtree
			if (now->left) {
				now = now->left;
			}
			else {//impossible condition
				cout << "interval doesn't exist" << endl;
				exit(-1);
			}
		}
		else {
			now->key.endT = newEndT;
			break;
		}
	}
}

//
//int main() {
//	
//	NewTree* a = new NewTree();
//	a->insert(Intv(1, 2, 1));
//	a->insert(Intv(3, 5, 1));
//	a->insert(Intv(6, 7, 1));
//	a->insert(Intv(8, 12, 1));
//	a->insert(Intv(13, 14, 1));
//	a->insert(Intv(15, 15, 1));
//	a->insert(Intv(16, 24, 1));
//	a->insert(Intv(25, 25, 1));
//	a->insert(Intv(26, 33, 1));
//	a->insert(Intv(33, 35, 1));
//	a->print();
//	cout << "==================================\n";
//	vector<Intv> result;
//	a->searchAndIntersect(result,2,23,3);
//	vector<Intv>::iterator iter= result.begin();
//	for (; iter != result.end(); iter++) {
//		cout << *iter << endl;
//	}
//	system("pause");
//}