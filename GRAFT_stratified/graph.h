#ifndef GRAPH_H_
#define GRAPH_H_
//without tuple look up
#include <fstream>
#include <exception>
#include <iostream>
#include <iterator>
#include <set>
#include <map>
#include "StringTokenizer.h"
#include "random.h"
#include "triangle.h"
#define BUCKETS 5
//#include "IQueue.h"

using namespace std;
using namespace boost;
typedef pair<int,int> EDGE;

struct EDGE_DEG {
	EDGE edge;
	unsigned int degree;
	unsigned int bucket;
};


struct prop_edge
{
	EDGE edge;
	bool validTC;
	int TC;//triangle count

	bool validDC;
	int DC;
};

struct prop_vertex
{
	int vertex;
	int adj_size;
	vector<int> adj;      //no order is maintained
};

struct compED{
	bool operator()(const EDGE_DEG& e1, const EDGE_DEG& e2)const{
		return e1.degree < e2.degree;
	}
}cmpED;

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
	map<int,double> symbolC;

	vector<prop_vertex> graphAdj;
	int vertexCount;

	map<int,int> vertexToIndex;//(vertex,index) pair
	bool vertexToIndexValid;

	vector<prop_edge> edge_list;
	int edgeCount;
	bool triCountDone;
	bool degreeCountDone;
	
	//bucket data structure
	vector<EDGE_DEG> edge_deg_list;
	vector<int> bucketBoundary;
	vector<int> bucketSize;
	vector<int> bucketSampleSize;
	//bucket data structure
	

	map<EDGE,int,cmpE> edgeToIndex;
	bool edgeToIndexValid;

	set<EDGE,cmpE> edgeSet;
	bool edgeSetValid;

	vector<int> edgeCW;
	bool edgeCWValid;

	//mcmc degree based
	int walkL;
	int randomizeFactor;
	int rFT;
	EDGE current;
	bool MCMCDInitiated;

	//mcmc triangle based
	bool MCMCTInitiated;

public:
	~graph_(){
		for(int i=0;i<graphAdj.size();i++){
			graphAdj[i].adj.clear();
		}
		graphAdj.clear();
		vertexToIndex.clear();
		edge_list.clear();
		edgeToIndex.clear();
		edgeSet.clear();
		edgeCW.clear();
		//cout<<"graph destructor\n";
	}

	graph_(graph_ *bg, const char* EorV, double reductionFactor){
		triCountDone=false;
		degreeCountDone=false;
    	edgeCWValid=false;
    	edgeToIndexValid=false;
    	edgeSetValid=false;
    	MCMCDInitiated=false;
    	MCMCTInitiated=false;
    	walkL=-1;
    	randomizeFactor=-1;
    	current=make_edge(0,0);
    	vertexToIndexValid=false;

    	if(strcmp(EorV,"edge")==0)
    	{
    		bg->load_reduced_edge_graph(reductionFactor,edge_list);
    		prepare_variable();
    	}else if(strcmp(EorV,"vertex")==0)
    	{
    		bg->load_reduced_vertex_graph(reductionFactor,edge_list);
    		prepare_variable();
    	}
    	else{
    		cout<<"error!!!!";
    		exit(1);
    	}
    	vertexToIndexValid=true;
    	//initiateEdgeToIndexMap();
    	initiateEdgeSet();
		//bucket data structure
		initiateBucketDataStructure();
	}

    graph_(const char* filename) {
    	triCountDone=false;
    	degreeCountDone=false;
    	edgeCWValid=false;
    	edgeToIndexValid=false;
    	MCMCDInitiated=false;
    	MCMCTInitiated=false;
    	walkL=-1;
    	randomizeFactor=-1;
    	current=make_edge(0,0);
    	vertexToIndexValid=false;

    	ifstream infile(filename, ios::in);
    	read_graph(infile);
    	vertexToIndexValid=true;
    	//initiateEdgeToIndexMap();
    	initiateEdgeSet();
		
		//bucket data structure
		initiateBucketDataStructure();
    }
	
	void initiateBucketDataStructure()
	{
		bucketBoundary.resize(BUCKETS+1);
		bucketSize.resize(BUCKETS+1);
		bucketSampleSize.resize(BUCKETS+1);
		unsigned int TotalED=0;
		unsigned int BucketSize=0;
        
        //Compute edge degree for each edge
		for (int i=0; i<edge_list.size(); i++) {
			EDGE_DEG e;
			e.edge=edge_list[i].edge;
			
			
			map<int,int>::const_iterator f,s;
			f=vertexToIndex.find(edge_list[i].edge.first);
			s=vertexToIndex.find(edge_list[i].edge.second);
			e.degree=graphAdj[f->second].adj.size();
			e.degree=e.degree + graphAdj[s->second].adj.size();
			TotalED=TotalED+e.degree; //get total edge degree
			
			edge_deg_list.push_back(e);
		}
		
		BucketSize=TotalED/BUCKETS; //each bucket has equal aggregated edge degree
		
		sort(edge_deg_list.begin(),edge_deg_list.end(),cmpED); //edges sorted on edge degree
		
		int t=0;
		int bn=1;
		bucketBoundary[0]=0;
		for (int i=0;i<edge_deg_list.size();i++)
		{
            //bucket is full and new edge has higher edge degree than the previous one
			if (t>=BucketSize and edge_deg_list[i].degree > edge_deg_list[i-1].degree) {
				bucketBoundary[bn]=i; //current bucket ends here
				bucketSize[bn]=bucketBoundary[bn]-bucketBoundary[bn-1];//adjact bucket size
				bn++;//we now fill the next bucket
				t=0;
			}
			edge_deg_list[i].bucket=bn;//assign an edge to a bucket
			t=t+edge_deg_list[i].degree;//total edge degree assigned to the current bucket
		}
		bucketBoundary[bn]=edge_deg_list.size();
		bucketSize[bn]=bucketBoundary[bn]-bucketBoundary[bn-1];
		printBucketDataStructure();
	}
	
	void printBucketDataStructure()
	{
		/*
         for (int i=0; i<edge_deg_list.size(); i++) {
			print_edge(edge_deg_list[i].edge);
			cout<<edge_deg_list[i].degree<<" "<<edge_deg_list[i].bucket <<"\n";
		}
        */
         
		cout<<"Bucket Boundary:: ";
		for (int i=0; i< bucketBoundary.size(); i++) {
			cout<<bucketBoundary[i]<<":";
		}
		cout<<"\n";
		
		cout<<"Bucket Sizes:: ";
		for (int i=1; i< bucketBoundary.size(); i++) {
			cout<<bucketSize[i]<<":";
		}
		cout<<"\n";
	}
	

    void load_reduced_vertex_graph(double reductionFactor,vector<prop_edge> &redge_list)
    {
    	double rand;
    	set<int> vertices;
    	set<int>::iterator it;
    	for(int i=0;i < graphAdj.size();i++){
    		if(graphAdj[i].adj_size==0) continue;
    		rand = random_uni01();
    		if(rand <= reductionFactor){
    			vertices.insert(graphAdj[i].vertex);
    		}
    	}

    	for(int i=0;i<edge_list.size();i++){
    		it=vertices.find(edge_list[i].edge.first);
    		if(it == vertices.end()) continue;
    		it=vertices.find(edge_list[i].edge.second);
    		if(it == vertices.end()) continue;

    		redge_list.push_back(edge_list[i]);
    	}


    }

    void load_reduced_edge_graph(double reductionFactor,vector<prop_edge> &redge_list){
    	double rand;
    	for(int i=0;i<edge_list.size();i++){
    		rand=random_uni01();
    		if(rand<=reductionFactor){
    			redge_list.push_back(edge_list[i]);
    			//redge_list[redge_list.size()-1].validTC=false;
    			//redge_list[redge_list.size()-1].validDC=false;
    		}
    	}
    }

    void prepare_variable()
    {
    	int largest_v=0;
    	for(int i=0;i<edge_list.size();i++){
    		edge_list[i].validDC=false;
    		edge_list[i].validTC=false;
    		if(edge_list[i].edge.second>largest_v) largest_v=edge_list[i].edge.second;
    	}
    	vertexCount = largest_v+1;
    	edgeCount=edge_list.size();
    	//make adjacency list
    	int vIndex=0;
    	map<int,int>::iterator it;
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
    	/*
    	graphAdj.resize(vertexCount);
    	for(int i=0;i<edge_list.size();i++)
    	{
       		graphAdj[edge_list[i].edge.first].adj.push_back(edge_list[i].edge.second);
       		graphAdj[edge_list[i].edge.second].adj.push_back(edge_list[i].edge.first);
    	}
    	for(int i=0;i<graphAdj.size();i++)
    	{
       		sort(graphAdj[i].adj.begin(),graphAdj[i].adj.end()); //pre-sort the adjacency list...needed for intersection operation
       		graphAdj[i].adj_size=graphAdj[i].adj.size();//pre-compute the adjacency list size
    		graphAdj[i].vertex=i;
    	}*/
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
    		e.validTC=false;
    		e.validDC=false;
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


    		//graphAdj[edge_list[i].edge.second].adj.push_back(edge_list[i].edge.first);
    	}
    	for(int i=0;i<graphAdj.size();i++)
    	{
    		sort(graphAdj[i].adj.begin(),graphAdj[i].adj.end()); //pre-sort the adjacency list...needed for intersection operation
    		graphAdj[i].adj_size=graphAdj[i].adj.size();//pre-compute the adjacency list size
    		//graphAdj[i].vertex=i;
    	}


    }

    /*int graph_triangle_count_exact()
    {
    	int tc=0;
    	for (int i=0;i<edge_list.size();i++){
    		tc=tc+triangle_count(edge_list[i].edge);
    	}
    	return tc;
    }*/

    int signature_count_app_total(int k, int sn, double p)
    {
    	int cumSignCount=0;
    	double rand;
    	if(k==4 and (sn == 2 or sn == 1 or sn == 5)){
    		for(int i=0;i<edge_list.size();i++)
    		{
    			rand=random_uni01();
    			if(rand <= p){
    				cumSignCount+=signature_count_edge(k,sn,edge_list[i].edge);
    			}
    		}
      		return cumSignCount/p;
    	}else if(k==4 and sn==4){
    		for(int i=0;i<edge_list.size();i++)
    		{
    			rand=random_uni01();
    			if(rand <= p){
      				cumSignCount+=signature_count_edge(k,sn,edge_list[i].edge);
    			}
    	   	}
     		return cumSignCount/(p*6);//for symmetry
     	}else if(k==4 and sn==3){
      		for(int i=0;i<edge_list.size();i++)
      		{
      			rand=random_uni01();
      			if(rand <= p){
      				cumSignCount+=signature_count_edge(k,sn,edge_list[i].edge);
      			}
      		}
    		return cumSignCount/(p*4);//for symmetry
    	}else if(k==4 and sn==6){
       		for(int i=0;i<edge_list.size();i++)
       		{
       			rand=random_uni01();
       			if(rand <= p){
       				cumSignCount+=signature_count_edge(k,sn,edge_list[i].edge);
       			}
       		}
       		return cumSignCount/(p*3);//for symmetry
    	}else if(k==3 and sn ==1){
    		return triMul3App(p)/3;
    	}
        return -1;
    }

    void print_stat(int edge_no){
    	//cout<<"edge#"<<edge_no<<": ";
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
		//cout<<symbolC[i]<<":";
		cout.flush();
    	}
       	cout<<"\n";
    }
    void signature_count_total_agg(int k){
    	for(int i=1;i<=29;i++){
       		symbolC[i]=0;
    	}

    	for(int i=1;i<=29;i++){
    		cout<<i<<":";
    	}
    	cout<<"\n";

    	for(int i=0;i<edge_list.size();i++)
    	{
    		//print_edge(edge_list[i].edge);
    		//cout<<"\n";
    		signature_count_edge_aggregated(k,edge_list[i].edge);
    		if ((i+1)%10000==0){
				print_stat(i+1);
			}
    	}



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



    	for(int i=1;i<=29;i++){
    		cout<<symbolC[i]<<":";
    	}
    	cout<<"\n";

    	/*for(int i=1 ;i<=29 ;i++){
    		cout<<i<<"\t"<<symbolC[i]<<"\n";
    	}*/
    }

	
    void signature_count_app_agg_bucket(int k, double p){
		vector<double> symbolCT;
		symbolCT.resize(30);
		for (int i=0; i<symbolCT.size(); i++) {
			symbolCT[i]=0;
		}
		//cout<<"CHECK1\n";
		set<int> sampleSet;
		
    	for (int i=1; i<bucketSize.size(); i++) {
			sampleSet.clear();
			bucketSampleSize[i]=bucketSize[i]*p;
			//cout << bucketSampleSize[i]<<" :samples\n";
			
			for(int j=1;j<=29;j++){
				symbolC[j]=0;
			}
			//cout<<"CHECK2\n";
			cout.flush();
			unsigned int rand;
			while(sampleSet.size()!=bucketSampleSize[i])
			{
				//cout<<"CHECK3\n";
				rand=boost_get_a_random_number(bucketBoundary[i-1],bucketBoundary[i]);
				if (sampleSet.find(rand)!=sampleSet.end()) {
					continue;
				}
				sampleSet.insert(rand);
				signature_count_edge_aggregated(k,edge_deg_list[rand].edge);
			}
			//if ((i+1)%10000==0){
			//cout.flush();
			
			//edge degree
			/*
			 map<int,int>::const_iterator f,s;
			 f=vertexToIndex.find(edge_list[i].edge.first);
			 s=vertexToIndex.find(edge_list[i].edge.second);
			 int ed=0;
			 ed=ed+graphAdj[f->second].adj.size();
			 ed=ed+graphAdj[s->second].adj.size();
			 cout<<ed<<":";*/
			//edge degree
			
			/*if ((i+1)%1000==0){
				print_stat(i+1);
				for (int i=1;i<=29;i++){
					symbolC[i]=0;
				}		
			}*/
			
    	
			//cout<<x<<" edges sampled\n";
		
		
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
		
			for(int j=1;j<=29;j++){
				symbolC[j]*=bucketSize[i];
				symbolC[j]/=bucketSampleSize[i];
			}
		
			for(int j=1;j<=29;j++){
				//cout<<symbolC[j]<<":";
				symbolCT[j]+=symbolC[j];
			}
			//cout<<"\n";
    	}
		for(int i=1;i<=29;i++){
			cout<<symbolCT[i]<<":";
		}
		cout<<"\n";
		
    }
	
	
    void signature_count_app_agg(int k, double p){
    	for(int i=1;i<=29;i++){
    		symbolC[i]=0;
    	}
    	double rand;
    	//int x=0;
    	for(int i=0;i<edge_list.size();i++)
    	{
    		rand=random_uni01();
    		if(rand <= p){
    			//x++;
    			signature_count_edge_aggregated(k,edge_list[i].edge);
    		}
			//if ((i+1)%10000==0){
				//cout.flush();
		
			//edge degree
			/*
			 map<int,int>::const_iterator f,s;
			 f=vertexToIndex.find(edge_list[i].edge.first);
			 s=vertexToIndex.find(edge_list[i].edge.second);
			 int ed=0;
			 ed=ed+graphAdj[f->second].adj.size();
			 ed=ed+graphAdj[s->second].adj.size();
			 cout<<ed<<":";*/
			//edge degree
		
			if ((i+1)%1000==0){
				//print_stat(i+1);
				/*for (int i=1;i<=29;i++){
					symbolC[i]=0;
				}*/		
			}
	
    	}
    	//cout<<x<<" edges sampled\n";


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


    	for(int i=1;i<=29;i++){
    		symbolC[i]/=p;
    	}

    	for(int i=1;i<=29;i++){
    		cout<<symbolC[i]<<":";
    	}
       	cout<<"\n";

    	/*for(int i=1 ;i<=29 ;i++){
    		symbolC[i]/=p;
    		cout<<i<<"\t"<<symbolC[i]<<"\n";
    	}*/
    }



    int signature_count_total(int k, int sn)
    {
    	int cumSignCount=0;
    	if(k==4 and (sn == 2 or sn == 1 or sn == 5)){
    		for(int i=0;i<edge_list.size();i++)
    		{
    			cumSignCount+=signature_count_edge(k,sn,edge_list[i].edge);
    		}
    		return cumSignCount;
    	}else if(k==4 and sn==4){
    		for(int i=0;i<edge_list.size();i++)
    		{
    			cumSignCount+=signature_count_edge(k,sn,edge_list[i].edge);
    		}
    		return cumSignCount/6;//for symmetry
    	}else if(k==4 and sn==3){
    		for(int i=0;i<edge_list.size();i++)
    		{
    			cumSignCount+=signature_count_edge(k,sn,edge_list[i].edge);
    		}
    		return cumSignCount/4;//for symmetry
    	}else if(k==4 and sn==6){
    		for(int i=0;i<edge_list.size();i++)
    		{
    			cumSignCount+=signature_count_edge(k,sn,edge_list[i].edge);
    		}
    		return cumSignCount/3;//for symmetry
    	}else if(k==3 and sn ==1){
    		return triMul3()/3;
    	}
        return -1;
    }

    bool findEdge(int a, int b){
    	return edgeSet.find(make_edge(a,b)) != edgeSet.end();
    }

    void signature_count_edge_aggregated(int k, EDGE e)
    {
    	//print_edge(e);
    	//cout<<"\n";
    	vector<int> v;
    	v.resize(k);

    	//01_,02_,03,04,12,13_,14,23,24_,34
    	//01_,02_,03,04,12,13_,14,23,24,34_
    	vector<bool> edgesT;
    	edgesT.resize(k*(k-1)/2);

    	//cout<<k*(k-1)/2<<"\n";
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
       	  	    				//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":11_1\n";
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
      	  	    				//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":10_1\n";
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
   		    					//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":11_2\n";
  	   			    			continue;
   		    				}

   		    			}
   					}
   				}
   				//no 11 end

   				//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<":10_2_undone\n";
   				for(int j=0;j<graphAdj[f->second].adj.size();j++){
    	    		v[3]=graphAdj[f->second].adj[j];
    	    		if (v[3]==s->first or v[2]==v[3]) continue;
    	    		//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]	<<":10_2_undone\n";

    	    		for(int k=j+1 ; k < graphAdj[f->second].adj.size() ; k++){
    	    			v[4]=graphAdj[f->second].adj[k];
    	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]	<< " "<<v[4]<<":10_2_undone\n";
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
      	  	    				//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":10_2\n";
      	  	    				continue;
      	  	    			}
      	  	    			symbolC[14]++;
      	  	    			continue;
      	  	    		}
      	  	    	//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":10_2_undone\n";
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
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":9\n";
    	  	    			continue;
    	  	    		}

    	  	    		if(edgesT[2]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[9] and edgesT[9]==false
    	  	    				and edgesT[4]!=false){
    	  	    			symbolC[12]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":12\n";
    	  	    			//print_edge(edgesT[4]->first);
    	  	    			//cout<<"\n";
    	  	    			continue;
    	  	    		}
    	  	    		//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":test\n";

    	  	    		if(edgesT[4]==edgesT[6]and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[9] and edgesT[9]==false
    	  	    				and edgesT[2]!=false and edgesT[3]!=false){
    	  	    			symbolC[18]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":18_1\n";
    	  	    			continue;
    	  	    		}else if(edgesT[4]==edgesT[9] and edgesT[9]==false
    	  	    				and edgesT[2]!=false and edgesT[3]!=false and edgesT[6]!=false and edgesT[7]!=false){
    	  	    			symbolC[27]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":18_1\n";
    	  	    			continue;
    	  	    		}else if(edgesT[4]==edgesT[9] and edgesT[9]==false
    	  	    				and edgesT[2]!=false and edgesT[3]!=false and
    	  	    				((edgesT[6]!=false and edgesT[7]==false)or (edgesT[6]==false and edgesT[7]!=false))){
    	  	    			symbolC[24]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":18_1\n";
    	  	    			continue;
    	  	    		}
    	  	    		else if(edgesT[4]==edgesT[6]and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[9] and edgesT[9]==false
    	  	    				and (edgesT[2]!=false or edgesT[3]!=false))
    	  	    		{
    	  	    			symbolC[13]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":13_1\n";
    	  	    			continue;
    	  	   			}


    	  	    		if(edgesT[4]==edgesT[2]and edgesT[2]==edgesT[3] and edgesT[3]==edgesT[9] and edgesT[9]==false
    	  	    		   		and edgesT[6]!=false and edgesT[7]!=false){
     	  	    			symbolC[20]++;
     	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":20_1\n";
    	  	    			continue;
    	  	    		}else if(edgesT[4]==edgesT[2]and edgesT[2]==edgesT[3] and edgesT[3]==edgesT[9] and edgesT[9]==false
    	  	    				and (edgesT[6]!=false or edgesT[7]!=false))
    	   	    		{
    	  	    			symbolC[16]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_1\n";
    	  	    			continue;
    	   	    		}

    	  	    		if(edgesT[4] != false and
    	   	    				edgesT[2]!=false and edgesT[3]!=false and
    	   	    				edgesT[6]!=false and edgesT[7]!=false)
    	   	    		{
    	  	    			if(edgesT[9]==false){
    	  	    				symbolC[28]++;
    	  	    				//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_1\n";
    	  	    				continue;
    	  	    			}
    	  	    			symbolC[29]++;
    	  	    			continue;
    	   	    		}else if(edgesT[9]==false and edgesT[4] != false and edgesT[6]!=false and edgesT[7]!=false and
    	   	    				((edgesT[2]==false and edgesT[3]!=false )
    	  	    						or (edgesT[3]==false and edgesT[2]!=false )))
    	   	    		{
    	  	    			symbolC[26]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_1\n";
    	  	    			continue;
    	   	    		}else if(edgesT[9]==false and edgesT[4] != false and
    	   	    				((edgesT[2]==false and edgesT[3]!=false and edgesT[6]!=false)
    	  	    						or (edgesT[3]==false and edgesT[2]!=false and edgesT[7]!=false)))
    	   	    		{
    	  	    			symbolC[23]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_1\n";
    	  	    			continue;
    	   	    		}else if(edgesT[9]==false and edgesT[4] != false and
    	   	    				edgesT[2]==false and edgesT[3]==false and
    	   	    				(edgesT[6]==false or edgesT[7]==false))
    	   	    		{
    	  	    			symbolC[17]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_1\n";
    	  	    			continue;
    	   	    		}else if(edgesT[9]==false and edgesT[4] == false and
    	   	    				((edgesT[7]!=false and edgesT[2]!=false and edgesT[6]==false and edgesT[3]==false)
    	   	    						or (edgesT[6]!=false and edgesT[3]!=false and edgesT[7]==false and edgesT[2]==false)))
    	   	    		{
    	  	    			symbolC[19]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_1\n";
    	  	    			continue;
    	   	    		}

    	  	    		if(edgesT[2]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[4] and edgesT[4]==false
    	  	    				and edgesT[9]!=false){
    	   	    			symbolC[15]++;
    	   	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":15\n";
    	  	    			//print_edge(edgesT[9]->first);
    	   	    			//cout<<"\n";
    	  	    			continue;
    	  	    		}
    	  	    		if(edgesT[2]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==false
    	  	    		   	and edgesT[9]!=false and edgesT[4]!=false){
    	  	    			symbolC[21]++;
    	   	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":15\n";
    	  	    			//print_edge(edgesT[9]->first);
    	   	    			//cout<<"\n";
       	  	    			continue;
    	  	    		}else if(edgesT[6]==edgesT[7] and edgesT[7]==false
        	  	    	   	and edgesT[9]!=false and edgesT[4]!=false and
        	  	    	   	((edgesT[2]!=false and edgesT[3]==false) or (edgesT[2]==false and edgesT[3]!=false))){
        	  	   			symbolC[25]++;
        	   	   			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":15\n";
        	  	   			//print_edge(edgesT[9]->first);
        	   	   			//cout<<"\n";
           	      			continue;
        	  	    	}
   			    	}


    	  	    	//getting adj list
    	  	    	//map<int,int>::const_iterator t;
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
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":9\n";
    	   	    			continue;
    	  	    		}

    	  	    		if(edgesT[4]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[8] and edgesT[8]==false
    	  	    				and edgesT[2]!=false){
    	  	      			symbolC[12]++;
    	  	      			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":12\n";
    	  	    			//print_edge(edgesT[2]->first);
    	  	    			//cout<<"\n";
    	   	    			continue;
    	  	    		}
    	  	    		//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":test2\n";

    	  	    		if(edgesT[2]==edgesT[3]and edgesT[3]==edgesT[7] and edgesT[7]==edgesT[8] and edgesT[8]==false
   	  	    				and edgesT[4]!=false and edgesT[6]!=false){
    	  	      			symbolC[18]++;
    	  	      			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":18_2\n";
    	   	    			continue;
    	  	    		}else if(edgesT[2]==edgesT[8] and edgesT[8]==false
       	  	    				and edgesT[4]!=false and edgesT[6]!=false and edgesT[3]!=false and edgesT[7]!=false){
        	  	      			symbolC[27]++;
        	  	      			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":18_2\n";
        	   	    			continue;
        	  	    	}else if(edgesT[2]==edgesT[8] and edgesT[8]==false
       	  	    				and edgesT[4]!=false and edgesT[6]!=false and
       	  	    				((edgesT[3]!=false and edgesT[7]==false) or (edgesT[3]==false and edgesT[7]!=false))){
        	  	      			symbolC[24]++;
        	  	      			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":18_2\n";
        	   	    			continue;
        	  	    	}
    	  	    		else if(edgesT[2]==edgesT[3]and edgesT[3]==edgesT[7] and edgesT[7]==edgesT[8] and edgesT[8]==false
   	  	    				and (edgesT[4]!=false or edgesT[6]!=false))
    	  	    		{
    	  	    			symbolC[13]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":13_2\n";
    	   	    			continue;
    	  	    		}


    	  	    		if(edgesT[2]==edgesT[4]and edgesT[4]==edgesT[6] and edgesT[6]==edgesT[8] and edgesT[8]==false
       	  	    				and edgesT[3]!=false and edgesT[7]!=false){
    	  	    			symbolC[20]++;
    	  	      			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":20_2\n";
    	  	    			continue;
    	  	    		}else if(edgesT[2]==edgesT[4]and edgesT[4]==edgesT[6] and edgesT[6]==edgesT[8] and edgesT[8]==false
    	  	    				and (edgesT[3]!=false or edgesT[7]!=false))
    	   	    		{
    	  	    			symbolC[16]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_2\n";
    	   	    			continue;
    	   	    		}


    	  	    		if(edgesT[2] != false and
    	  	    				edgesT[4]!=false and edgesT[6]!=false and
    	  	    				edgesT[3]!=false and edgesT[7]!=false)
    	   	    		{
    	  	    			if(edgesT[8]==false){
    	  	    				symbolC[28]++;
    	  	    				//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_1\n";
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
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_1\n";
    	  	    			continue;
    	   	    		}else if(edgesT[8]==false and edgesT[2] != false and
    	   	    				((edgesT[4]==false and edgesT[6]!=false and edgesT[3]!=false)
    	  	    						or (edgesT[6]==false and edgesT[4]!=false and edgesT[7]!=false)))
    	   	    		{
    	  	    			symbolC[23]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_1\n";
    	  	    			continue;
    	   	    		}else if(edgesT[8]==false and edgesT[2] != false and
    	   	    				edgesT[4]==false and edgesT[6]==false and
    	   	    				(edgesT[3]==false or edgesT[7]==false))
    	   	    		{
    	  	    			symbolC[17]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_2\n";
    	   	    			continue;
    	   	    		}else if(edgesT[2] == false and edgesT[8]==false and
    	   	    				((edgesT[7]!=false and edgesT[4]!=false and edgesT[3]==false and edgesT[6]==false) or
    	   	    						(edgesT[3]!=false and edgesT[6]!=false and edgesT[7]==false and edgesT[4]==false)))
    	   	    		{
    	  	    			symbolC[19]++;
    	  	    			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":16_2\n";
    	   	    			continue;
    	   	    		}


    	  	    		if(edgesT[4]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==edgesT[2] and edgesT[2]==false
    	  	    				and edgesT[8]!=false){
    	   	      			symbolC[15]++;
    	   	      			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":15\n";
    	   	      			//print_edge(edgesT[8]->first);
    	  	    			//cout<<"\n";
    	   	      			continue;
    	  	    		}
    	  	    		if(edgesT[4]==edgesT[3]and edgesT[3]==edgesT[6] and edgesT[6]==edgesT[7] and edgesT[7]==false
    	  	    		   		and edgesT[8]!=false and edgesT[2]!=false){
    	   	      			symbolC[21]++;
    	   	      			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":15\n";
    	   	      			//print_edge(edgesT[8]->first);
    	   	      			//cout<<"\n";
    	   	      			continue;
    	  	    		}else if(edgesT[3]==edgesT[7] and edgesT[7]==false
    	  	    		   		and edgesT[8]!=false and edgesT[2]!=false and
    	  	    		   		((edgesT[4]!=false and edgesT[6]==false) or (edgesT[4]==false and edgesT[6]!=false))){
    	   	      			symbolC[25]++;
    	   	      			//cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<":15\n";
    	   	      			//print_edge(edgesT[8]->first);
    	   	      			//cout<<"\n";
    	   	      			continue;
    	  	    		}
    	  	    	}
    	      	}
    	    }
       	}

    }

    int signature_count_edge(int k, int sn, EDGE e)
    {
       	if(k==4){
    		if(sn==1){
    			map<int,int>::const_iterator f,s;
    			f=vertexToIndex.find(e.first);
       			s=vertexToIndex.find(e.second);
       			if(f == vertexToIndex.end() or s == vertexToIndex.end())
       			{
       				cout<<"the vertex does not exist\n";
    	    		exit(1);
       			}
       			if(find(graphAdj[f->second].adj.begin(), graphAdj[f->second].adj.end(), e.second)==graphAdj[f->second].adj.end()
       					or find(graphAdj[s->second].adj.begin(), graphAdj[s->second].adj.end(), e.first)==graphAdj[s->second].adj.end()){
       				print_edge(e);
    				cout<<"the edge does not exist\n";
    				exit(1);
    			}
       			vector<int> v;
       			v.resize(graphAdj[f->second].adj.size()+graphAdj[s->second].adj.size());
       			vector<int>::iterator it,itt;
   		    	it=set_intersection(graphAdj[f->second].adj.begin(),graphAdj[f->second].adj.end(),graphAdj[s->second].adj.begin(),graphAdj[s->second].adj.end(),v.begin());
    			itt=v.begin();
    			int sn1=0;
    			while(itt!=it){
    				int a1 = (*itt);
    				map<int,int>::const_iterator a1a=vertexToIndex.find(a1);
    				for(int j=0;j<graphAdj[a1a->second].adj.size();j++)
    				{
    				    int a2=graphAdj[a1a->second].adj[j];
    				    if (a2==s->first or a2==f->first) continue;
    				    if (edgeToIndex.find(make_edge(a2,s->first)) != edgeToIndex.end()) continue;
    				    if (edgeToIndex.find(make_edge(a2,f->first)) != edgeToIndex.end()) continue;
    				    sn1++;
    				}
    				itt++;
    			}
   		    	return sn1;
    		}else if(sn==3){
    			map<int,int>::const_iterator f,s;
    			f=vertexToIndex.find(e.first);
    			s=vertexToIndex.find(e.second);
    			if(f == vertexToIndex.end() or s == vertexToIndex.end())
    			{
    	    		cout<<"the vertex does not exist\n";
   		    		exit(1);
    			}
    			if(find(graphAdj[f->second].adj.begin(), graphAdj[f->second].adj.end(), e.second)==graphAdj[f->second].adj.end()
    			   		or find(graphAdj[s->second].adj.begin(), graphAdj[s->second].adj.end(), e.first)==graphAdj[s->second].adj.end()){
    				print_edge(e);
    				cout<<"the edge does not exist\n";
    	    		//return 0; //the edge does not exist
    				exit(1);
    			}

    			int sn3=0;
    			for (int i=0;i<graphAdj[f->second].adj.size();i++){
    			    int a1=graphAdj[f->second].adj[i];
    			    if (a1==s->first) continue;

    			    for(int j=0;j<graphAdj[s->second].adj.size();j++){
    			    	int a2=graphAdj[s->second].adj[j];
    			    	if (a2==f->first) continue;
    			    	if (edgeToIndex.find(make_edge(a2,a1)) != edgeToIndex.end()){
    			    		if (edgeToIndex.find(make_edge(a2,f->first)) != edgeToIndex.end()) continue;
    			    		if (edgeToIndex.find(make_edge(s->first,a1)) != edgeToIndex.end()) continue;
    			    		sn3++;
    			    	}
    			    }
    			}
   		    	return sn3;
    		}else if(sn==2){
    			map<int,int>::const_iterator f,s;
    			f=vertexToIndex.find(e.first);
    			s=vertexToIndex.find(e.second);
    			if(f == vertexToIndex.end() or s == vertexToIndex.end())
    			{
    	    		cout<<"the vertex does not exist\n";
   		    		exit(1);
    			}
    			if(find(graphAdj[f->second].adj.begin(), graphAdj[f->second].adj.end(), e.second)==graphAdj[f->second].adj.end()
    			   		or find(graphAdj[s->second].adj.begin(), graphAdj[s->second].adj.end(), e.first)==graphAdj[s->second].adj.end()){
    				print_edge(e);
    				cout<<"the edge does not exist\n";
    	    		//return 0; //the edge does not exist
    				exit(1);
    			}
    			vector<int> v;
    			v.resize(graphAdj[f->second].adj.size()+graphAdj[s->second].adj.size());
    	    	vector<int>::iterator it,itt1,itt2;
   		    	it=set_intersection(graphAdj[f->second].adj.begin(),graphAdj[f->second].adj.end(),graphAdj[s->second].adj.begin(),graphAdj[s->second].adj.end(),v.begin());
   		    	itt1=v.begin();
   		    	int sn2=0;
   		    	while(itt1!=it){
   		    		map<int,int>::const_iterator a1a=vertexToIndex.find((*itt1));
   		    		itt2=itt1;
   		    		itt2++;
   		    		while(itt2!=it){
   		    			map<int,int>::const_iterator a2a=vertexToIndex.find((*itt2));
   		    			///you can chose the option not to use this code(output will be same)
   		    			if(edgeToIndex.find(make_edge(a2a->first,a1a->first))==edgeToIndex.end()){
   		    			   	sn2++;
   		    			}
   		      			itt2++;
   		    		}
   		    		itt1++;
   		    	}
   		    	return sn2;
    		}else if(sn==4){
    			map<int,int>::const_iterator f,s;
    			f=vertexToIndex.find(e.first);
       			s=vertexToIndex.find(e.second);
       			if(f == vertexToIndex.end() or s == vertexToIndex.end())
       			{
       				cout<<"the vertex does not exist\n";
    	    		exit(1);
       			}
       			if(find(graphAdj[f->second].adj.begin(), graphAdj[f->second].adj.end(), e.second)==graphAdj[f->second].adj.end()
       					or find(graphAdj[s->second].adj.begin(), graphAdj[s->second].adj.end(), e.first)==graphAdj[s->second].adj.end()){
       				print_edge(e);
    				cout<<"the edge does not exist\n";
    				exit(1);
    			}
       			vector<int> v;
       			v.resize(graphAdj[f->second].adj.size()+graphAdj[s->second].adj.size());
       			vector<int>::iterator it,itt1,itt2;
   		    	it=set_intersection(graphAdj[f->second].adj.begin(),graphAdj[f->second].adj.end(),graphAdj[s->second].adj.begin(),graphAdj[s->second].adj.end(),v.begin());
   		    	itt1=v.begin();
   		    	int sn4=0;
   		    	while(itt1!=it){
   		    		map<int,int>::const_iterator a1a=vertexToIndex.find((*itt1));
   		    		itt2=itt1;
   		    		itt2++;
   		    		while(itt2!=it){
   		    			map<int,int>::const_iterator a2a=vertexToIndex.find((*itt2));
   		    			//print_edge(make_edge(a2a->first,a1a->first));
   		    			if(edgeToIndex.find(make_edge(a2a->first,a1a->first)) != edgeToIndex.end()){
   		    				sn4++;
   		    			}
   		    			itt2++;
   		    		}
   		    		itt1++;
   		    	}
   		    	return sn4;
    		}
    		else if(sn==6)
    		{
    			map<int,int>::const_iterator f,s;
    			f=vertexToIndex.find(e.first);
      			s=vertexToIndex.find(e.second);
      			if(f == vertexToIndex.end() or s == vertexToIndex.end())
       			{
      				cout<<"the vertex does not exist\n";
    	    		exit(1);
       			}
      			if(find(graphAdj[f->second].adj.begin(), graphAdj[f->second].adj.end(), e.second)==graphAdj[f->second].adj.end()
      					or find(graphAdj[s->second].adj.begin(), graphAdj[s->second].adj.end(), e.first)==graphAdj[s->second].adj.end()){
      				print_edge(e);
    				cout<<"the edge does not exist\n";
    				exit(1);
    			}
      			int sn6=0;
      			for (int i=0;i<graphAdj[f->second].adj.size();i++){
      				int a1=graphAdj[f->second].adj[i];
      			    if (a1==s->first) continue;
      			    for(int j=i+1;j<graphAdj[f->second].adj.size();j++){
      			    	int a2=graphAdj[f->second].adj[j];
      			    	if (a2==s->first) continue;
      			    	if (edgeToIndex.find(make_edge(a2,a1)) != edgeToIndex.end()) continue;
      			    	if (edgeToIndex.find(make_edge(a1,s->first)) != edgeToIndex.end()) continue;
      			    	if (edgeToIndex.find(make_edge(a2,s->first)) != edgeToIndex.end()) continue;
      			    	sn6++;
      			    }
      			}

      			for (int i=0;i<graphAdj[s->second].adj.size();i++){
      				int a1=graphAdj[s->second].adj[i];
       			    if (a1==f->first) continue;
       			    for(int j=i+1;j<graphAdj[s->second].adj.size();j++){
       			    	int a2=graphAdj[s->second].adj[j];
      			    	if (a2==f->first) continue;
      			    	if (edgeToIndex.find(make_edge(a2,a1)) != edgeToIndex.end()) continue;
      			    	if (edgeToIndex.find(make_edge(a1,f->first)) != edgeToIndex.end()) continue;
      			    	if (edgeToIndex.find(make_edge(a2,f->first)) != edgeToIndex.end()) continue;
      			    	sn6++;
       			    }
      			}

      			return sn6;
    		}else if(sn==5)
    		{
    			map<int,int>::const_iterator f,s;
     			f=vertexToIndex.find(e.first);
     			s=vertexToIndex.find(e.second);
       			if(f == vertexToIndex.end() or s == vertexToIndex.end())
       			{
       	    		cout<<"the vertex does not exist\n";
       	    		exit(1);
    			}
       			if(find(graphAdj[f->second].adj.begin(), graphAdj[f->second].adj.end(), e.second)==graphAdj[f->second].adj.end()
       					or find(graphAdj[s->second].adj.begin(), graphAdj[s->second].adj.end(), e.first)==graphAdj[s->second].adj.end()){
       				print_edge(e);
       				cout<<"the edge does not exist\n";
    				exit(1);
       			}
       			int sn5=0;

       			for (int i=0;i<graphAdj[f->second].adj.size();i++){
       			    int a1=graphAdj[f->second].adj[i];
       			    if (a1==s->first) continue;

       			    for(int j=0;j<graphAdj[s->second].adj.size();j++){
      			    	int a2=graphAdj[s->second].adj[j];
      			    	if (a2==f->first or a1==a2) continue;
      			    	if (edgeToIndex.find(make_edge(a2,a1)) != edgeToIndex.end()) continue;
      			    	if (edgeToIndex.find(make_edge(a1,s->first)) != edgeToIndex.end()) continue;
      			    	if (edgeToIndex.find(make_edge(a2,f->first)) != edgeToIndex.end()) continue;
      			    	sn5++;
       			    }
       			}
       			return sn5;
       		}

    	}
        return -1;
    }

    int triangle_count(EDGE e)
    {
    	map<int,int>::const_iterator f,s;
    	f=vertexToIndex.find(e.first);
       	s=vertexToIndex.find(e.second);
    	if(f == vertexToIndex.end() or s == vertexToIndex.end())
    	{
    		cout<<"the vertex does not exist\n";
    		//return -1;//the vertex does not exist
    		exit(1);
    	}
    	if(find(graphAdj[f->second].adj.begin(), graphAdj[f->second].adj.end(), e.second)==graphAdj[f->second].adj.end()
    			or find(graphAdj[s->second].adj.begin(), graphAdj[s->second].adj.end(), e.first)==graphAdj[s->second].adj.end()){
    		print_edge(e);
    		cout<<"the edge does not exist\n";
    		//return 0; //the edge does not exist
    		exit(1);
    	}
    	vector<int> v;
    	v.resize(graphAdj[f->second].adj.size()+graphAdj[s->second].adj.size());
    	vector<int>::iterator it;
    	it=set_intersection(graphAdj[f->second].adj.begin(),graphAdj[f->second].adj.end(),graphAdj[s->second].adj.begin(),graphAdj[s->second].adj.end(),v.begin());
    	return it-v.begin();
    }

    int degree_count(EDGE e)
    {
    	map<int,int>::const_iterator f,s;
    	f=vertexToIndex.find(e.first);
       	s=vertexToIndex.find(e.second);
       	if(f == vertexToIndex.end() or s == vertexToIndex.end())
      	{
       		cout<<"the vertex does not exist\n";
      		exit(1);
      	}
       	if(find(graphAdj[f->second].adj.begin(), graphAdj[f->second].adj.end(), e.second)==graphAdj[f->second].adj.end()
       			or find(graphAdj[s->second].adj.begin(), graphAdj[s->second].adj.end(), e.first)==graphAdj[s->second].adj.end()){
    	   	print_edge(e);
     		cout<<"the edge does not exist\n";
     		exit(1);
    	}
       	return graphAdj[f->second].adj.size()+graphAdj[s->second].adj.size();
    }


    void precompute_all_DC()
    {
    	if(degreeCountDone==true) return;
    	for(int i=0;i<edge_list.size();i++)
    	{
       		if(edge_list[i].validDC==false){
       			edge_list[i].DC=degree_count(edge_list[i].edge);
    	   		edge_list[i].validDC=true;
     		}
    	}
    	degreeCountDone=true;
    }

    void precompute_all_TC()
    {
    	if(triCountDone==true) return;

    	for(int i=0;i<edge_list.size();i++)
		{
    		if(edge_list[i].validTC==false){
    			edge_list[i].TC=triangle_count(edge_list[i].edge);
    			edge_list[i].validTC=true;
    		}
		}
    	triCountDone=true;
    	/*for(int i=0;i<edge_list.size();i++)
    	{
    		if(edge_list[i].validTC==false)
    		{
    			cout<<"error!!! TC not counted\n";
    		}
    		cout<<"edge: ("<<edge_list[i].edge.first<<","<<edge_list[i].edge.second<<") TC: "<<edge_list[i].TC<<'\n';
    	}*/
    }

    int uniformEdgeSampleTC()//returns triangle count of a uniformly sampled edge
    {
    	int rand = boost_get_a_random_number(0,edge_list.size());
    	//print_edge(edge_list[rand].edge);
    	//cout<< triangle_count(edge_list[rand].edge)<<endl;
    	return triangle_count(edge_list[rand].edge);
    }

    EDGE uniformEdgeSample()
    {
    	int rand = boost_get_a_random_number(0,edge_list.size());
    	return edge_list[rand].edge;
    }

    triangle* SampleTriangle(const char* method)
    {
    	EDGE e;
    	if(strcmp(method,"UE")==0){
    		e=uniformEdgeSample();
    		while (triangle_count(e)==0){
    			e=uniformEdgeSample();
    		}
    	}else if(strcmp(method,"WE")==0){
    		e=weightedEdgeSample();
    		while (triangle_count(e)==0){
    			e=weightedEdgeSample();
    		}
    	}else if(strcmp(method,"MCMCD")==0){
    		//cout<<"1";
    		initiateMCMCDegree(5,6);
    		//cout<<"2";
    		MCMCDegreeEdgeSample();
    		//cout<<"3";
    		while (triangle_count(current)==0){
    			MCMCDegreeEdgeSample();
    			//cout<<"4";
    		}
    		e=current;
    	}else if (strcmp(method,"MCMCT")==0){
    		initiateMCMCTriangle(5,6);
    		MCMCTriangleEdgeSample();
    		while (triangle_count(current)==0){
    			MCMCTriangleEdgeSample();
    		}
    		e=current;
    	}
    	else{
    		return NULL;
    	}
    	//cout<<e.first<<e.second;
    	return sampleTriangleOfEdge(e);
    }

    triangle* sampleTriangleOfEdge(EDGE e)
    {
    	map<int,int>::const_iterator f,s;
    	f=vertexToIndex.find(e.first);
       	s=vertexToIndex.find(e.second);
       	if(f == vertexToIndex.end() or s == vertexToIndex.end())
      	{
       		cout<<"the vertex does not exist\n";
     		exit(1);
      	}
       	if(find(graphAdj[f->second].adj.begin(), graphAdj[f->second].adj.end(), e.second)==graphAdj[f->second].adj.end()
       			or find(graphAdj[s->second].adj.begin(), graphAdj[s->second].adj.end(), e.first)==graphAdj[s->second].adj.end()){
    	   	print_edge(e);
      		cout<<"the edge does not exist\n";
      		exit(1);
    	}
    	vector<int> v;
       	v.resize(graphAdj[f->second].adj.size()+graphAdj[s->second].adj.size());
       	vector<int>::iterator it;
       	it=set_intersection(graphAdj[f->second].adj.begin(),graphAdj[f->second].adj.end(),graphAdj[s->second].adj.begin(),graphAdj[s->second].adj.end(),v.begin());
       	int rand=boost_get_a_random_number(0,it-v.begin());
       	triangle *t=new triangle(v[rand],e.first,e.second);
       	return t;
    }


    void initiateWeightedSampleTC()
    {
    	if(edgeCWValid==true) return;

    	precompute_all_TC();
       	edgeCW.resize(edge_list.size());
       	edgeCW[0]=edge_list[0].TC;
       	for(int i=1;i<edge_list.size();i++)
       	{
       		edgeCW[i]=edgeCW[i-1]+edge_list[i].TC;
       	}

       	edgeCWValid=true;
    }

    void initiateEdgeSet()
    {
    	if(edgeSetValid==true) return;
    	for(int i=0;i<edge_list.size();i++){
      		edgeSet.insert(edge_list[i].edge);
    	}
    	edgeSetValid=true;
    }

    void initiateEdgeToIndexMap()
    {
    	if(edgeToIndexValid==true) return;
    	for(int i=0;i<edge_list.size();i++){
    		edgeToIndex[edge_list[i].edge]=i;
    	}
    	cout<<"in initiateEdgeToIndexMap";
    	edgeToIndexValid=true;

    	/*for(int i=0;i<edge_list.size();i++)
    	{
    		print_edge(edge_list[edgeToIndex[edge_list[i].edge]].edge);
    		print_edge(edge_list[i].edge);
    		cout<<" "<<i<<"\n";
    	}*/
    }

    void initiateMCMCDegree(int walkL_,int randomizeFactor_)
    {
    	if(MCMCDInitiated==true){
    		return;
    	}
    	walkL=walkL_;
    	randomizeFactor=randomizeFactor_;
    	rFT=randomizeFactor_;
    	precompute_all_DC();
    	initiateEdgeToIndexMap();
    	current=uniformEdgeSample();
    	MCMCDInitiated=true;
    }

    void initiateMCMCTriangle(int walkL_,int randomizeFactor_)
    {
    	if(MCMCTInitiated==true){
    	    return;
    	}
    	walkL=walkL_;
       	randomizeFactor=randomizeFactor_;
       	rFT=randomizeFactor_;
       	precompute_all_DC();
       	initiateEdgeToIndexMap();
       	current=uniformEdgeSample();
       	MCMCTInitiated=true;
    }


    void MCMCTriangleEdgeSample()
    {
       	if(rFT==0){
       		rFT=randomizeFactor;
       		current=uniformEdgeSample();
       	}
      	else{
       		rFT--;
       	}
       	//EDGE next;
       	for(int i=0;i<walkL;i++){
       		//cout<<"5";
        	selectNextTB();
       	}
    }

    void MCMCDegreeEdgeSample()
    {
    	if(rFT==0){
    		rFT=randomizeFactor;
    		current=uniformEdgeSample();
    	}
    	else{
    		rFT--;
    	}
    	//EDGE next;
    	for(int i=0;i<walkL;i++){
    		//cout<<"5";
    		selectNextDB();
    	}
    }


    void selectNextTB(){
        EDGE next;
        int ei=edgeToIndex[current];
        int cdc=edge_list[ei].DC;
        int rand=boost_get_a_random_number(0,cdc);
        if(rand < graphAdj[vertexToIndex[edge_list[ei].edge.first]].adj_size){
        	next= make_edge(edge_list[ei].edge.first, graphAdj[vertexToIndex[edge_list[ei].edge.first]].adj[rand]);
       	}else{
       		rand=rand-graphAdj[vertexToIndex[edge_list[ei].edge.first]].adj_size;
       		next= make_edge(edge_list[ei].edge.second, graphAdj[vertexToIndex[edge_list[ei].edge.second]].adj[rand]);
       	}

        if(edge_list[ei].validTC==false){
        	edge_list[ei].TC=triangle_count(edge_list[ei].edge);
        	//print_edge(edge_list[ei].edge);
          	edge_list[ei].validTC=true;
        }
        int ctc=edge_list[ei].TC;

        ei=edgeToIndex[next];
        if(edge_list[ei].validTC==false){
            edge_list[ei].TC=triangle_count(edge_list[ei].edge);
            //print_edge(edge_list[ei].edge);
           	edge_list[ei].validTC=true;
        }
        int ndc=edge_list[ei].TC;
        rand=boost_get_a_random_number(0,cdc);
       	if(rand<ndc){
       	    current=next;
       	}
    }


    void selectNextDB(){
    	EDGE next;
    	int ei=edgeToIndex[current];
    	//print_edge(current);
    	int cdc=edge_list[ei].DC;
    	//cout<<cdc;
    	int rand=boost_get_a_random_number(0,cdc);
    	//cout << " "<<rand<<"\n";
    	if(rand < graphAdj[vertexToIndex[edge_list[ei].edge.first]].adj_size){
    		next= make_edge(edge_list[ei].edge.first, graphAdj[vertexToIndex[edge_list[ei].edge.first]].adj[rand]);
    	}else{
    		rand=rand-graphAdj[vertexToIndex[edge_list[ei].edge.first]].adj_size;
    		next= make_edge(edge_list[ei].edge.second, graphAdj[vertexToIndex[edge_list[ei].edge.second]].adj[rand]);
    	}
    	//print_edge(next);
    	//cdc=edge_list[ei].DC;
    	int ndc=edge_list[edgeToIndex[next]].DC;
    	/*if(cdc <= ndc){
    		current=next;
    		//return current;
    	}else{
    		rand=boost_get_a_random_number(0,cdc);
    		if(rand<ndc){
    			current=next;
    			//return current;
    		}
    	}//return current;*/
    	rand=boost_get_a_random_number(0,cdc);
    	if(rand<ndc){
    	    current=next;
    	}
    }

    int weightedEdgeSampleTC()//returns triangle count of a weighted(based on triangle count) sampled edge
    {
    	initiateWeightedSampleTC();
    	int rand = boost_get_a_random_number(1,edgeCW[edgeCW.size()-1]+1);
    	//binary search
    	int x=0,y=edgeCW.size()-1;
    	int mid;
    	int edgeI=-1;

    	while(y!=x+1){
    		mid=(x+y)/2;
    		if(edgeCW[mid]<rand){
    			x=mid;
    		}else if(edgeCW[mid]>=rand and edgeCW[mid-1]<rand){
    			edgeI=mid;
    			break;
    		}else{
    			y=mid;
    		}
    	}
    	if (edgeI==-1)
    	{
    		if(edgeCW[y]>=rand and edgeCW[y-1]<rand){
    			edgeI=y;
    		}else{
    			edgeI=x;
    		}
    	}
    	//print_edge(edge_list[edgeI].edge);
    	//cout<< triangle_count(edge_list[edgeI].edge)<<endl;
    	return triangle_count(edge_list[edgeI].edge);
    }

    EDGE weightedEdgeSample()
    {
    	initiateWeightedSampleTC();
    	int rand = boost_get_a_random_number(1,edgeCW[edgeCW.size()-1]+1);
    	//binary search
    	int x=0,y=edgeCW.size()-1;
    	int mid;
    	int edgeI=-1;

    	while(y!=x+1){
    		mid=(x+y)/2;
    		if(edgeCW[mid]<rand){
    			x=mid;
    	   	}else if(edgeCW[mid]>=rand and edgeCW[mid-1]<rand){
     			edgeI=mid;
     			break;
    	   	}else{
    			y=mid;
    	   	}
    	}
      	if (edgeI==-1)
      	{
       		if(edgeCW[y]>=rand and edgeCW[y-1]<rand){
       			edgeI=y;
    	   	}else{
    			edgeI=x;
    	   	}
      	}
      	//print_edge(edge_list[edgeI].edge);
      	//cout<< triangle_count(edge_list[edgeI].edge)<<endl;
       	return edge_list[edgeI].edge;
    }

    bool isTriPresent(triangle *t)
    {
    	vector<int> v;
    	t->getVertices(v);
    	map<EDGE,int,cmpE>::const_iterator cit;
    	initiateEdgeToIndexMap();
    	cit=edgeToIndex.find(make_edge(v[0],v[1]));
    	if(cit==edgeToIndex.end()) return false;
    	cit=edgeToIndex.find(make_edge(v[0],v[2]));
    	if(cit==edgeToIndex.end()) return false;
    	cit=edgeToIndex.find(make_edge(v[1],v[2]));
    	if(cit==edgeToIndex.end()) return false;
    	return true;
    }

    bool isTriPresent(const vector<int> &v)
    {
    	map<EDGE,int,cmpE>::const_iterator cit;
       	initiateEdgeToIndexMap();
       	cit=edgeToIndex.find(make_edge(v[0],v[1]));
       	if(cit==edgeToIndex.end()) return false;
       	cit=edgeToIndex.find(make_edge(v[0],v[2]));
       	if(cit==edgeToIndex.end()) return false;
       	cit=edgeToIndex.find(make_edge(v[1],v[2]));
       	if(cit==edgeToIndex.end()) return false;
       	return true;
    }

    void print_adjacency_list()
    {
      	cout<<"-------graph begin-------\n";
      	for(int i=0;i<graphAdj.size();i++)
      	{
       		cout<<"vertex: "<<graphAdj[i].vertex<<" degree: "<<graphAdj[i].adj_size<<"---- ";
       		for(int j=0;j<graphAdj[i].adj.size();j++)
       		{
       			cout<<graphAdj[i].adj[j]<<", ";
       		}
       		cout<<"\n";
      	}
       	cout<<"-------graph end-------\n";
    }

    EDGE make_edge(int a, int b)
    {
    	if (a<b) return EDGE(a,b);
     	return EDGE(b,a);
    }

    void print_edge(EDGE e)
    {
    	cout<<"("<<e.first<<","<<e.second<<") ";
    }

    long triMul3()
    {
    	initiateWeightedSampleTC();
    	return edgeCW[edgeCW.size()-1];
    }

    long triMul3App(double p){
    	double rand;
       	long TCMul3=0;
       	for(int i=0;i<edge_list.size();i++)
   		{
       		rand = random_uni01();
        	if(rand <= p)
       		{
       			if(edge_list[i].validTC==false){
       				edge_list[i].TC=triangle_count(edge_list[i].edge);
       				edge_list[i].validTC=true;
   	    		}
       			TCMul3+=edge_list[i].TC;
        	}
       	}
       	return TCMul3/p;
    }

    int getEdgeCount(){return edgeCount;}

    int getVertexCount(){return vertexToIndex.size();}

};
#endif
