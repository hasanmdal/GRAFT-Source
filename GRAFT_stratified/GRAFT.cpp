#include <string.h>
#include <vector>
#include "graph.h"
#include "time_tracker.h"
using namespace std;

char* datafile;
double sample_fraction;

void print_usage(char *prog) {
  cerr<<"Usage: "<<prog<<" -d data-file -sf sample_fraction(0...1)"<<endl;
  exit(0);
}

void parse_args(int argc, char* argv[]) {
  if(argc<5) {
    print_usage(argv[0]);
  }
  for (int i=1; i < argc; ++i){
    if (strcmp(argv[i], "-d") == 0){
    	datafile=argv[++i];
    }
    else if(strcmp(argv[i], "-sf") == 0){
    	sample_fraction=atof(argv[++i]);
   }
  }
}//end parse_args()

int main(int argc, char *argv[]) {
	typedef struct induction ind;
	parse_args(argc, argv);
	graph_* g=new graph_(datafile);
	cout << "Data File: "<<datafile<<"\n";
    cout<<"Sample Fraction: "<<sample_fraction<<"\n";

	cout<<"Total Edges: "<<g->getEdgeCount()<<"\n";
	cout.flush();
	cout<<"Total Vertices: "<<g->getVertexCount()<<"\n";
	cout.flush();
	double fs=sample_fraction;

    //cout << "exact: ";
	//g->signature_count_app_agg(5,1);
	cout<<"Simple::\n";
    g->signature_count_app_agg(5,fs);
	cout<<"Strafied with five buckets::\n";
    g->signature_count_app_agg_bucket(5,fs);
    delete g;
}
