#ifndef GRAPH_H_
#define GRAPH_H_
#include <fstream>
#include <exception>
#include <iostream>
#include <iterator>
#include <set>
#include <map>
#include "StringTokenizer.h"
#include "random.h"
	
using namespace std;
using namespace boost;
typedef pair<int,int> EDGE;

struct prop_edge //Legacy: this structure was suppose to be more complex, but letter those extra properties were not used.
{
	EDGE edge;
};

struct prop_vertex
{
	int vertex;
	int adj_size;
	vector<int> adj;      //no order is maintained
};


struct cmpE{
	bool operator()(const EDGE& e1,const EDGE& e2)const{
		if(e1.first<e2.first){
			return true;
		}else if(e1.first>e2.first){
			return false;
		}else{
			if(e1.second<e2.second){
				return true;
			}
			return false;
		}
	}
};
class graph_{
public:
	typedef pair<int,int> EDGE;
private:
	map<int,double> symbolC;//counts of graphlet are stored here

	vector<prop_vertex> graphAdj;//graph adjacency list and edgelist both are used to make our operations efficient.
	int vertexCount;

	map<int,int> vertexToIndex;//(vertex,index) pair
	bool vertexToIndexValid;

	
	//Improvement possible: edge_list and edgeSet can be merged, one of the structures is enough for our porpose
	vector<prop_edge> edge_list;
	int edgeCount;
	
	set<EDGE,cmpE> edgeSet;
	bool edgeSetValid;

public:
	~graph_(){
		for(int i=0;i<graphAdj.size();i++){
			graphAdj[i].adj.clear();
		}
		graphAdj.clear();
		vertexToIndex.clear();
		edge_list.clear();
		edgeSet.clear();
	}


    graph_(const char* filename) {
    	vertexToIndexValid=false;
		edgeSetValid=false;
		 
    	ifstream infile(filename, ios::in);
    	read_graph(infile);
    	vertexToIndexValid=true;
    	initiateEdgeSet();
    }

	
    void read_graph(ifstream& datafile)
    {

    	/*
    	 * v1\tv2\n
    	 * Graph should be undirected. No duplicated edge must be present.
    	 * For best result, vertex labels should start from 0 and end to n.
    	 * All numbers between 0 to n should represent a valid node of graph.
    	 * File should end with a newline '\n'
    	 */
    	std::string oneline;
    	int largest_v=0;

    	while(1)								//read all the edges in a vector
    	{
    		std::getline(datafile,oneline);
    		if (oneline.length() < 1) break;
    		StringTokenizer strtok = StringTokenizer(oneline,"\t");
    		int v1 = strtok.nextIntToken();
    		int v2 = strtok.nextIntToken();
    		if (v1>largest_v) largest_v=v1;
    		if (v2>largest_v) largest_v=v2;
    		prop_edge e;
    		e.edge= make_edge(v1,v2);
    		edge_list.push_back(e);
    	}
    	vertexCount = largest_v+1;
    	edgeCount = edge_list.size();
    	//make adjacency list
    	int vIndex=0;
    	map<int,int>::iterator it;
    	//graphAdj.resize(vertexCount);
    	for(int i=0;i<edge_list.size();i++)
    	{
    		it=vertexToIndex.find(edge_list[i].edge.first);
    		if(it == vertexToIndex.end()){
    			prop_vertex a;
    			graphAdj.push_back(a);
    			graphAdj[vIndex].vertex=edge_list[i].edge.first;
    			vertexToIndex.insert(pair<int,int>(edge_list[i].edge.first,vIndex));
    			vIndex++;
    		}
    		int ti=vertexToIndex.find(edge_list[i].edge.first)->second;
    		if(edge_list[i].edge.first!=graphAdj[ti].vertex) {
    			cout<<"error read_graph!!!\n";
    			exit(1);
    		}
    		graphAdj[ti].adj.push_back(edge_list[i].edge.second);

    		it=vertexToIndex.find(edge_list[i].edge.second);
    		if(it == vertexToIndex.end()){
       			prop_vertex a;
       			graphAdj.push_back(a);
       			graphAdj[vIndex].vertex=edge_list[i].edge.second;
    		    vertexToIndex.insert(pair<int,int>(edge_list[i].edge.second,vIndex));
    		    vIndex++;
    		}
    		ti=vertexToIndex.find(edge_list[i].edge.second)->second;
    		if(edge_list[i].edge.second!=graphAdj[ti].vertex) {
    			cout<<"error read_graph!!!\n";
       			exit(1);
    		}
    		graphAdj[ti].adj.push_back(edge_list[i].edge.first);
    	}
    	for(int i=0;i<graphAdj.size();i++)
    	{
    		sort(graphAdj[i].adj.begin(),graphAdj[i].adj.end()); //pre-sort the adjacency list...needed for intersection operation
    		graphAdj[i].adj_size=graphAdj[i].adj.size();//pre-compute the adjacency list size
    	}


    }
    
    void print_stat(int edge_no){
    	cout<<"edge#"<<edge_no<<": ";
	cout.flush();
    	vector<float> t;
    	t.resize(30);

    	t[1]=symbolC[1]/2;
    	t[2]=symbolC[2]/6;
    	t[3]=symbolC[3];
    	t[4]=symbolC[4]/3;
    	t[5]=symbolC[5]/4;
    	t[6]=symbolC[6]/2;
    	t[7]=symbolC[7]/4;
    	t[8]=symbolC[8]/12;
    	t[9]=symbolC[9]/2;
    	t[10]=symbolC[10];
    	t[12]=symbolC[12]/2;
    	t[13]=symbolC[13]/4;
    	t[18]=symbolC[18]/8;
    	t[15]=symbolC[15]/10;
    	t[16]=symbolC[16]/4;
    	t[20]=symbolC[20]/12;
    	t[21]=symbolC[21]/2;
    	t[23]=symbolC[23]/12;
    	t[17]=symbolC[17]/4;
    	t[28]=symbolC[28]/12;
    	t[26]=symbolC[26]/8;
    	t[19]=symbolC[19]/4;
    	t[29]=symbolC[29]/120;
    	t[14]=symbolC[14]/2;
    	t[11]=symbolC[11]/4;
    	t[22]=symbolC[22]/12;
    	t[27]=symbolC[27]/8;
    	t[24]=symbolC[24]/4;
    	t[25]=symbolC[25]/8;

    	for(int i=1;i<=29;i++){
    		cout<<t[i]<<":";
		cout.flush();
    	}
       	cout<<"\n";
    }
    
    void signature_count_app_agg(double p){//counts all the graphlets
    	for(int i=1;i<=29;i++){//initiating counter
    		symbolC[i]=0;
    	}
    	double rand;
    	for(int i=0;i<edge_list.size();i++)
    	{
    		rand=random_uni01();//generate uniform random number
    		if(rand <= p){
    			signature_count_edge_aggregated(edge_list[i].edge);
    		}
			//used if intermediate results are needed
			/*if ((i+1)%1==0){  
				print_stat(i+1);
			}*/

    	}
    	

		//normalizing the grphlet counts (due to automorphism)
    	symbolC[1]/=2;
    	symbolC[2]/=6;
    	symbolC[4]/=3;
    	symbolC[5]/=4;
    	symbolC[6]/=2;
    	symbolC[7]/=4;
    	symbolC[8]/=12;
    	symbolC[9]/=2;
    	symbolC[12]/=2;
    	symbolC[13]/=4;
    	symbolC[18]/=8;
    	symbolC[15]/=10;
    	symbolC[16]/=4;
    	symbolC[20]/=12;
    	symbolC[21]/=2;
    	symbolC[23]/=12;
    	symbolC[17]/=4;
    	symbolC[28]/=12;
    	symbolC[26]/=8;
    	symbolC[19]/=4;
    	symbolC[29]/=120;
    	symbolC[14]/=2;
    	symbolC[11]/=4;
    	symbolC[22]/=12;
    	symbolC[27]/=8;
    	symbolC[24]/=4;
    	symbolC[25]/=8;

		//normalizing the grphlet counts (due to sampling)
    	for(int i=1;i<=29;i++){
    		symbolC[i]/=p;
    	}
		
		//output the approximated count
    	for(int i=1;i<=29;i++){
    		cout<<symbolC[i]<<":";
    	}
       	cout<<"\n";
    }


    bool findEdge(int a, int b){
    	return edgeSet.find(make_edge(a,b)) != edgeSet.end();
    }

    void signature_count_edge_aggregated(EDGE e)//counts all the graphlets embedded with edge "e" as First aligned edge (FAE) and vertex count form 3 to 5
    {
		int k=5;//number of vertices of largest graphlet we considure
    	
		vector<int> v;
    	v.resize(k);

    	//01_,02_,03,04,12,13_,14,23,24_,34
    	//01_,02_,03,04,12,13_,14,23,24,34_
    	vector<bool> edgesT;
    	edgesT.resize(k*(k-1)/2);

       	if(k==5){
       		map<int,int>::const_iterator f,s;
    	    f=vertexToIndex.find(e.first);
   			s=vertexToIndex.find(e.second);
   			if(f == vertexToIndex.end() or s == vertexToIndex.end())
   			{
   				cout<<"the vertex does not exist\n";
        		exit(1);
   			}
   			v[0]=f->first;//v[0] contains the first vertex of edge
   			v[1]=s->first;//v[1] contains the second vertex of edge
   			if(find(graphAdj[f->second].adj.begin(), graphAdj[f->second].adj.end(), e.second)==graphAdj[f->second].adj.end()
   					or find(graphAdj[s->second].adj.begin(), graphAdj[s->second].adj.end(), e.first)==graphAdj[s->second].adj.end()){
    	  		print_edge(e);
    	  		cout<<"the edge does not exist\n";
    			exit(1);
   			}
   			//new
   			for (int i=0;i<graphAdj[f->second].adj.size();i++){
   				v[2]=graphAdj[f->second].adj[i];
   				if (v[2]==s->first) continue;

   				//3_nodes
   				if(findEdge(v[1],v[2])==false){
   					symbolC[1]++;
   				}else{
   					symbolC[2]++;
   				}

   				//end_3_nodes

   				//no 11
   				for(int j=i+1; j<graphAdj[f->second].adj.size();j++){
   					v[3]=graphAdj[f->second].adj[j];
   					if(v[3]==s->first) continue;

   					//4_node_star
   					if(findEdge(v[1],v[2]) == false
   							and findEdge(v[1],v[3]) == false
   							and findEdge(v[2],v[3]) == false){
   						symbolC[4]++;
   					}
   					//end_4_node_star

   					for(int k=j+1; k<graphAdj[f->second].adj.size();k++){
   						v[4]=graphAdj[f->second].adj[k];
   						if(v[4]==s->first) continue;
   						edgesT[4]=findEdge(v[1],v[2]);
   						edgesT[5]=findEdge(v[1],v[3]);
       	    			edgesT[6]=findEdge(v[1],v[4]);
       	    			edgesT[7]=findEdge(v[2],v[3]);
       	  	    		edgesT[8]=findEdge(v[2],v[4]);
       	  	    		edgesT[9]=findEdge(v[3],v[4]);

       	  	    		//can be more simplified
       	  	    		if (edgesT[7] ==false and edgesT[8]==false){
       	  	    			if(edgesT[4]==false and edgesT[5]==false and
           	  	    				edgesT[6]==false and edgesT[9]==false){
       	  	    				symbolC[11]++;
       	  	        			continue;
       	  	    			}
       	  	    		}
   					}
   				}
   				//no 11 end

   				for(int j=0;j<graphAdj[s->second].adj.size();j++){
    	    		v[3]=graphAdj[s->second].adj[j];
    	    		if (v[3]==f->first or v[2]==v[3]) continue;

    	    		//4_node_except_4
    	    		edgesT[2]=findEdge(v[0],v[3]);
    	    		edgesT[3]=findEdge(v[1],v[2]);
    	    		edgesT[5]=findEdge(v[2],v[3]);

    	    		if (edgesT[2]==false and edgesT[3]==false and edgesT[5] == false){
    	    			symbolC[3]++;
    	    		}else if(edgesT[2]==false and edgesT[3]== false and edgesT[5] != false){
    	    			symbolC[5]++;
    	    		}else if(edgesT[5] == false and
    	    				((edgesT[2] == false and edgesT[3] != false) or
    	    						(edgesT[2] != false and edgesT[3] == false)) ){
    	    			symbolC[6]++;
    	    		}else if(edgesT[5] != false and
    	    				((edgesT[2] == false and edgesT[3] != false) or
    	    						(edgesT[2] != false and edgesT[3] == false)) ){
    	    			symbolC[7]++;
    	    		}else if(edgesT[2] != false and
    	    				edgesT[3] != false and edgesT[5] != false){
    	    			symbolC[8]++;
    	    		}
    	    		//end_4_node_except_4

    	    		for(int k=j+1;k<graphAdj[s->second].adj.size();k++){
    	    			v[4]=graphAdj[s->second].adj[k];
    	    			if (v[4]==v[0] or v[4]==v[1] or v[4]==v[3]) continue;
    	    			edgesT[2]=findEdge(v[0],v[3]);
    	    			edgesT[3]=findEdge(v[0],v[4]);
    	    			edgesT[4]=findEdge(v[1],v[2]);
    	    			edgesT[7]=findEdge(v[2],v[3]);
      	  	    		edgesT[8]=findEdge(v[2],v[4]);
      	  	    		edgesT[9]=findEdge(v[3],v[4]);
      	  	    		if (edgesT[2]==edgesT[3] and edgesT[3]==edgesT[7]and
      	  	    				edgesT[7] ==edgesT[8] and
      	  	    				edgesT[8]==edgesT[9] and edgesT[9]==false){
      	  	    			if(edgesT[4] == false){
      	  	    				symbolC[10]++;
      	  	    				continue;
      	  	    			}
      	  	    			symbolC[14]++;
      	  	    			continue;
      	  	    		}
    	    		}
    	    	}
   			}

   			for (int i=0;i<graphAdj[s->second].adj.size();i++){
   				v[2]=graphAdj[s->second].adj[i];
   				if (v[2]==f->first) continue;

   				//3_nodes
   				if(findEdge(v[0],v[2])==false){
   					symbolC[1]++;
   				}else{
   					symbolC[2]++;
   				}
   				//end_3_nodes


   				//no 11
   				for(int j=i+1; j<graphAdj[s->second].adj.size();j++){
   					v[3]=graphAdj[s->second].adj[j];
   					if(v[3]==f->first) continue;

   					//4_node_star
   					if(findEdge(v[0],v[2]) == false
   							and findEdge(v[0],v[3]) == false
   							and findEdge(v[2],v[3]) == false){
   						symbolC[4]++;
   					}
   					//end_4_node_star

   					for(int k=j+1; k<graphAdj[s->second].adj.size();k++){
   						v[4]=graphAdj[s->second].adj[k];
   						if(v[4]==f->first) continue;
   						edgesT[1]=findEdge(v[0],v[2]);
   						edgesT[2]=findEdge(v[0],v[3]);
   						edgesT[3]=findEdge(v[0],v[4]);
   		    			edgesT[7]=findEdge(v[2],v[3]);
   		    			edgesT[8]=findEdge(v[2],v[4]);
   		    			edgesT[9]=findEdge(v[3],v[4]);

   		    			//can be more simplified
   		    			if (edgesT[7] ==edgesT[8] and edgesT[8]==false){
   		    				if(edgesT[1]==edgesT[2] and edgesT[2]==edgesT[3] and
   		    						edgesT[3]==edgesT[9]and edgesT[9]==false){
   		    					symbolC[11]++;
   		    	    			continue;
   		    				}

   		    			}
   					}
   				}
   				//no 11 end

   				for(int j=0;j<graphAdj[f->second].adj.size();j++){
    	    		v[3]=graphAdj[f->second].adj[j];
    	    		if (v[3]==s->first or v[2]==v[3]) continue;
    	    		
    	    		for(int k=j+1 ; k < graphAdj[f->second].adj.size() ; k++){
    	    			v[4]=graphAdj[f->second].adj[k];
    	    			if (v[4]==v[0] or v[4]==v[1] or v[4]==v[3]) continue;
    	    			edgesT[1]=findEdge(v[0],v[2]);
    	    			edgesT[5]=findEdge(v[1],v[3]);
    	    			edgesT[6]=findEdge(v[1],v[4]);
    	    			edgesT[7]=findEdge(v[2],v[3]);
      	  	    		edgesT[8]=findEdge(v[2],v[4]);
      	  	    		edgesT[9]=findEdge(v[3],v[4]);
      	  	    		if (edgesT[5]==edgesT[6] and edgesT[6]==edgesT[7]and
      	  	    				edgesT[7] ==edgesT[8] and
      	  	    				edgesT[8]==edgesT[9] and edgesT[9]==false){
      	  	    			if(edgesT[1] == false){
      	  	    				symbolC[10]++;
      	  	    				continue;
      	  	    			}
      	  	    			symbolC[14]++;
      	  	    			continue;
      	  	    		}
      	  	    	}
    	    	}
   			}


   			//End_new



    	    for (int i=0;i<graphAdj[f->second].adj.size();i++){
    	    	v[2]=graphAdj[f->second].adj[i];
    	    	if (v[2]==s->first) continue;
    		    for(int j=0;j<graphAdj[s->second].adj.size();j++){
    		    	v[3]=graphAdj[s->second].adj[j];
    	  	    	if (v[3]==f->first or v[2]==v[3]) continue;
    	  	    	//getting adj list
    	  	    	map<int,int>::const_iterator t;
    	  	    	t=vertexToIndex.find(v[2]);
        	    	if(t == vertexToIndex.end())
        	    	{
        	    		cout<<"the vertex does not exist\n";
        	    		exit(1);
    	  	    	}
    	  	    	for(int k=0;k<graphAdj[t->second].adj.size();k++){
    	  	    		v[4]=graphAdj[t->second].adj[k];
    	  	    		if (v[4]==v[0] or v[4]==v[1] or v[4]==v[3]) continue;
    	  	    		//01_,02_,03,04,12,13_,14,23,24_,34
    	  	    		edgesT[2]=findEdge(v[0],v[3]);
    	  	    		edgesT[3]=findEdge(v[0],v[4]);
    	  	    		edgesT[4]=findEdge(v[1],v[2]);
    	  	    		edgesT[6]=findEdge(v[1],v[4]);
    	  	    		edgesT[7]=findEdge(v[2],v[3]);
    	  	    		edgesT[9]=findEdge(v[3],v[4]);

    	  	    		if (edgesT[2]==edgesT[3] and edgesT[3]==edgesT[9] and edgesT[9]==false and
    	  	    				edgesT[4]!=false and edgesT[6]!=false and edgesT[7]!=false){
    	  	    			symbolC[22]++;
    	  	    			continue;
    	  	    		}

    	  	    		if (edgesT[2]==edgesT[3] and edgesT[3]==edgesT[4]and edgesT[4]==edgesT[6]and edgesT[6] ==edgesT[7] and edgesT[7]==edgesT[9] and edgesT[9]==false){
    	  	    			symbolC[9]++;
    	  	    			continue;
    	  	    		}

    	  	    		if(edgesT[2]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[9] and edgesT[9]==false
    	  	    				and edgesT[4]!=false){
    	  	    			symbolC[12]++;
    	  	    			continue;
    	  	    		}
    	  	    	
    	  	    		if(edgesT[4]==edgesT[6]and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[9] and edgesT[9]==false
    	  	    				and edgesT[2]!=false and edgesT[3]!=false){
    	  	    			symbolC[18]++;
    	  	    			continue;
    	  	    		}else if(edgesT[4]==edgesT[9] and edgesT[9]==false
    	  	    				and edgesT[2]!=false and edgesT[3]!=false and edgesT[6]!=false and edgesT[7]!=false){
    	  	    			symbolC[27]++;
    	  	    			continue;
    	  	    		}else if(edgesT[4]==edgesT[9] and edgesT[9]==false
    	  	    				and edgesT[2]!=false and edgesT[3]!=false and
    	  	    				((edgesT[6]!=false and edgesT[7]==false)or (edgesT[6]==false and edgesT[7]!=false))){
    	  	    			symbolC[24]++;
    	  	    			continue;
    	  	    		}
    	  	    		else if(edgesT[4]==edgesT[6]and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[9] and edgesT[9]==false
    	  	    				and (edgesT[2]!=false or edgesT[3]!=false))
    	  	    		{
    	  	    			symbolC[13]++;
    	  	    			continue;
    	  	   			}


    	  	    		if(edgesT[4]==edgesT[2]and edgesT[2]==edgesT[3] and edgesT[3]==edgesT[9] and edgesT[9]==false
    	  	    		   		and edgesT[6]!=false and edgesT[7]!=false){
     	  	    			symbolC[20]++;
     	  	    			continue;
    	  	    		}else if(edgesT[4]==edgesT[2]and edgesT[2]==edgesT[3] and edgesT[3]==edgesT[9] and edgesT[9]==false
    	  	    				and (edgesT[6]!=false or edgesT[7]!=false))
    	   	    		{
    	  	    			symbolC[16]++;
    	  	    			continue;
    	   	    		}

    	  	    		if(edgesT[4] != false and
    	   	    				edgesT[2]!=false and edgesT[3]!=false and
    	   	    				edgesT[6]!=false and edgesT[7]!=false)
    	   	    		{
    	  	    			if(edgesT[9]==false){
    	  	    				symbolC[28]++;
    	  	    				continue;
    	  	    			}
    	  	    			symbolC[29]++;
    	  	    			continue;
    	   	    		}else if(edgesT[9]==false and edgesT[4] != false and edgesT[6]!=false and edgesT[7]!=false and
    	   	    				((edgesT[2]==false and edgesT[3]!=false )
    	  	    						or (edgesT[3]==false and edgesT[2]!=false )))
    	   	    		{
    	  	    			symbolC[26]++;
    	  	    			continue;
    	   	    		}else if(edgesT[9]==false and edgesT[4] != false and
    	   	    				((edgesT[2]==false and edgesT[3]!=false and edgesT[6]!=false)
    	  	    						or (edgesT[3]==false and edgesT[2]!=false and edgesT[7]!=false)))
    	   	    		{
    	  	    			symbolC[23]++;
    	  	    			continue;
    	   	    		}else if(edgesT[9]==false and edgesT[4] != false and
    	   	    				edgesT[2]==false and edgesT[3]==false and
    	   	    				(edgesT[6]==false or edgesT[7]==false))
    	   	    		{
    	  	    			symbolC[17]++;
    	  	    			continue;
    	   	    		}else if(edgesT[9]==false and edgesT[4] == false and
    	   	    				((edgesT[7]!=false and edgesT[2]!=false and edgesT[6]==false and edgesT[3]==false)
    	   	    						or (edgesT[6]!=false and edgesT[3]!=false and edgesT[7]==false and edgesT[2]==false)))
    	   	    		{
    	  	    			symbolC[19]++;
    	  	    			continue;
    	   	    		}

    	  	    		if(edgesT[2]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[4] and edgesT[4]==false
    	  	    				and edgesT[9]!=false){
    	   	    			symbolC[15]++;
    	   	    			continue;
    	  	    		}
    	  	    		if(edgesT[2]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==false
    	  	    		   	and edgesT[9]!=false and edgesT[4]!=false){
    	  	    			symbolC[21]++;
    	   	    			continue;
    	  	    		}else if(edgesT[6]==edgesT[7] and edgesT[7]==false
        	  	    	   	and edgesT[9]!=false and edgesT[4]!=false and
        	  	    	   	((edgesT[2]!=false and edgesT[3]==false) or (edgesT[2]==false and edgesT[3]!=false))){
        	  	   			symbolC[25]++;
        	   	   			continue;
        	  	    	}
   			    	}


    	  	    	t=vertexToIndex.find(v[3]);
        	    	if(t == vertexToIndex.end())
        	    	{
        	    		cout<<"the vertex does not exist\n";
        	    		exit(1);
    	  	    	}
    	  	    	for(int k=0;k<graphAdj[t->second].adj.size();k++){
    	  	    		v[4]=graphAdj[t->second].adj[k];
    	  	    		if (v[4]==v[0] or v[4]==v[1] or v[4]==v[2]) continue;
    	  	    		//01_,02_,03,04,12,13_,14,23,24,34_
    	  	    		edgesT[2]=findEdge(v[0],v[3]);
    	  	    		edgesT[3]=findEdge(v[0],v[4]);
    	  	    		edgesT[4]=findEdge(v[1],v[2]);
    	  	    		edgesT[6]=findEdge(v[1],v[4]);
    	  	    		edgesT[7]=findEdge(v[2],v[3]);
    	  	    		edgesT[8]=findEdge(v[2],v[4]);

    	  	    		if (edgesT[4]==edgesT[6] and edgesT[6]==edgesT[8] and edgesT[8]==false and
    	  	    				edgesT[2]!=false and edgesT[3]!=false and edgesT[7]!=false){
    	  	    		   	 symbolC[22]++;
    	  	    		   	 continue;
    	  	    		}

    	  	    		if (edgesT[2]==edgesT[3] and edgesT[3]==edgesT[4]and edgesT[4]==edgesT[6]and edgesT[6] ==edgesT[7] and edgesT[7]==edgesT[8] and edgesT[8]==false){
    	  	    			symbolC[9]++;
    	  	    			continue;
    	  	    		}

    	  	    		if(edgesT[4]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[8] and edgesT[8]==false
    	  	    				and edgesT[2]!=false){
    	  	      			symbolC[12]++;
    	  	      			continue;
    	  	    		}
    	  	    		
    	  	    		if(edgesT[2]==edgesT[3]and edgesT[3]==edgesT[7] and edgesT[7]==edgesT[8] and edgesT[8]==false
   	  	    				and edgesT[4]!=false and edgesT[6]!=false){
    	  	      			symbolC[18]++;
    	  	      			continue;
    	  	    		}else if(edgesT[2]==edgesT[8] and edgesT[8]==false
       	  	    				and edgesT[4]!=false and edgesT[6]!=false and edgesT[3]!=false and edgesT[7]!=false){
        	  	      			symbolC[27]++;
        	  	      			continue;
        	  	    	}else if(edgesT[2]==edgesT[8] and edgesT[8]==false
       	  	    				and edgesT[4]!=false and edgesT[6]!=false and
       	  	    				((edgesT[3]!=false and edgesT[7]==false) or (edgesT[3]==false and edgesT[7]!=false))){
        	  	      			symbolC[24]++;
        	  	      			continue;
        	  	    	}
    	  	    		else if(edgesT[2]==edgesT[3]and edgesT[3]==edgesT[7] and edgesT[7]==edgesT[8] and edgesT[8]==false
   	  	    				and (edgesT[4]!=false or edgesT[6]!=false))
    	  	    		{
    	  	    			symbolC[13]++;
    	  	    			continue;
    	  	    		}


    	  	    		if(edgesT[2]==edgesT[4]and edgesT[4]==edgesT[6] and edgesT[6]==edgesT[8] and edgesT[8]==false
       	  	    				and edgesT[3]!=false and edgesT[7]!=false){
    	  	    			symbolC[20]++;
    	  	      			continue;
    	  	    		}else if(edgesT[2]==edgesT[4]and edgesT[4]==edgesT[6] and edgesT[6]==edgesT[8] and edgesT[8]==false
    	  	    				and (edgesT[3]!=false or edgesT[7]!=false))
    	   	    		{
    	  	    			symbolC[16]++;
    	  	    			continue;
    	   	    		}


    	  	    		if(edgesT[2] != false and
    	  	    				edgesT[4]!=false and edgesT[6]!=false and
    	  	    				edgesT[3]!=false and edgesT[7]!=false)
    	   	    		{
    	  	    			if(edgesT[8]==false){
    	  	    				symbolC[28]++;
    	  	    				continue;
    	  	    			}
    	  	    			symbolC[29]++;
    	  	    			continue;
    	   	    		}else if(edgesT[8]==false and edgesT[2] != false
    	   	    				and edgesT[3]!=false and edgesT[7]!=false and
    	   	    				((edgesT[4]==false and edgesT[6]!=false)
    	  	    						or (edgesT[6]==false and edgesT[4]!=false)))
    	   	    		{
    	  	    			symbolC[26]++;
    	  	    			continue;
    	   	    		}else if(edgesT[8]==false and edgesT[2] != false and
    	   	    				((edgesT[4]==false and edgesT[6]!=false and edgesT[3]!=false)
    	  	    						or (edgesT[6]==false and edgesT[4]!=false and edgesT[7]!=false)))
    	   	    		{
    	  	    			symbolC[23]++;
    	  	    			continue;
    	   	    		}else if(edgesT[8]==false and edgesT[2] != false and
    	   	    				edgesT[4]==false and edgesT[6]==false and
    	   	    				(edgesT[3]==false or edgesT[7]==false))
    	   	    		{
    	  	    			symbolC[17]++;
    	  	    			continue;
    	   	    		}else if(edgesT[2] == false and edgesT[8]==false and
    	   	    				((edgesT[7]!=false and edgesT[4]!=false and edgesT[3]==false and edgesT[6]==false) or
    	   	    						(edgesT[3]!=false and edgesT[6]!=false and edgesT[7]==false and edgesT[4]==false)))
    	   	    		{
    	  	    			symbolC[19]++;
    	  	    			continue;
    	   	    		}


    	  	    		if(edgesT[4]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[2] and edgesT[2]==false
    	  	    				and edgesT[8]!=false){
    	   	      			symbolC[15]++;
    	   	      			continue;
    	  	    		}
    	  	    		if(edgesT[4]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==false
    	  	    		   		and edgesT[8]!=false and edgesT[2]!=false){
    	   	      			symbolC[21]++;
    	   	      			continue;
    	  	    		}else if(edgesT[3]==edgesT[7] and edgesT[7]==false
    	  	    		   		and edgesT[8]!=false and edgesT[2]!=false and
    	  	    		   		((edgesT[4]!=false and edgesT[6]==false) or (edgesT[4]==false and edgesT[6]!=false))){
    	   	      			symbolC[25]++;
    	   	      			continue;
    	  	    		}
    	  	    	}
    	      	}
    	    }
       	}

    }


    EDGE make_edge(int a, int b)
    {
    	if (a<b) return EDGE(a,b);
     	return EDGE(b,a);
    }

    void initiateEdgeSet()
    {
    	if(edgeSetValid==true) return;
    	for(int i=0;i<edge_list.size();i++){
      		edgeSet.insert(edge_list[i].edge);
    	}
    	edgeSetValid=true;
    }

    void print_edge(EDGE e)
    {
    	cout<<"("<<e.first<<","<<e.second<<") ";
    }

    int getEdgeCount(){return edgeCount;}

    int getVertexCount(){return vertexToIndex.size();}

};
#endif
