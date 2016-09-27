/*
 * triangle.h
 *
 *  Created on: Nov 23, 2011
 *      Author: mmrahman
 */

#ifndef TUPLE_SET_H_
#define TUPLE_SET_H_

#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

struct node{
	int vertex;
	vector<node>* next;
};

typedef struct node node;
typedef vector<node>::iterator nodeScaner;

bool myComp (node i,node j) { return (i.vertex < j.vertex); }

class tuple_set
{
private:
	unsigned int tuple_length;
	node my_set;

public:
	tuple_set(int tl)
	{
		tuple_length=tl;
		my_set.vertex=-1;
		my_set.next=new vector<node>;
		cout<<"In constructor\n";
	}

	bool insert (vector<int>& tup_){
		bool prev_existing=true;

		if(tup_.size() < tuple_length){
			cout<<"tuple size invalid\n";
			exit(1);
		}
		vector<int> tup;
		for(int i=0;i<tuple_length;i++){
			tup.push_back(tup_[i]);
			//cout<<tup_[i]<<" ";
		}
		//cout<<"\n";

		sort(tup.begin(),tup.end());

		node* tnode;
		tnode = &my_set;
		nodeScaner tnext;

		for(int i=0;i<tup.size();i++){
			tnext=find_(tnode->next->begin(),tnode->next->end(),tnode->next->end(),tup[i]);
			if(tnext==tnode->next->end()){
				prev_existing=false;
				node t;
				t.vertex=tup[i];
				t.next=new vector<node>;
				tnode->next->push_back(t);
				sort(tnode->next->begin(),tnode->next->end(),myComp); //keeping the vector of nodes in sorted order. so that we can perform binary search
				tnode=&*(find_(tnode->next->begin(),tnode->next->end(),tnode->next->end(),tup[i]));
				continue;
			}
			tnode=&*tnext;
			//cout<<tup[i]<<"\n";
		}
		//cout<<"test5\n";
		tup.clear();
		return !prev_existing;
	}

	nodeScaner find_(nodeScaner b, nodeScaner e, nodeScaner end, int val){
		/*nodeScaner i=b;
		while(i!=e){
			if(i->vertex==val) return i;
			i++;
		}
		return e;*/
		if(b==e) return end;
		int t=(e-b)/2;
		nodeScaner i = b+t;
		if(i->vertex==val) return i;
		else if (i->vertex > val) return find_(b,i,end,val);
		else return find_(i+1,e,end,val);
		//binary find
	}

	bool find(vector<int>& tup_){
		if(tup_.size() < tuple_length){
			cout<<"tuple size invalid\n";
			exit(1);
		}

		vector<int> tup;
		for(int i=0;i<tuple_length;i++){
			tup.push_back(tup_[i]);
			//cout<<tup_[i]<<" ";
		}
		//cout<<"\n";

		sort(tup.begin(),tup.end());

		/*if(find(tup)==true){
			return false;
		}*/
		node* tnode;
		tnode = &my_set;
		nodeScaner tnext;
		//cout<<tuple_length<<" test5_2\n";
		for(int i=0;i<tup.size();i++){
			//cout<<"bi-search1\n";
			tnext=find_(tnode->next->begin(),tnode->next->end(),tnode->next->end(),tup[i]);
			//cout<<"bi-search2\n";
			if(tnext==tnode->next->end()){
				//cout<<"test6\n";
				tup.clear();
				return false;
			}
			tnode=&*tnext;
			//cout<<tup[i]<<"\n";
		}
		//cout<<"test7\n";
		tup.clear();
		return true;


	}
	~tuple_set(){
	}
};

#endif /* TUPLE_SET_H_ */
