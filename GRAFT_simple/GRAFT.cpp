#include <string.h>
#include <vector>
#include "graph.h"
#include "time_tracker.h"
using namespace std;

char* datafile;
double sample_fraction;

void print_usage(char *prog) {
  cerr<<"Usage: "<<prog<<" -d data-file -sf sample_fraction(0...1)"<<endl;
  cerr<<"input file format: \n\tEach line will contain an edge \n\tEdges are undirected, unweighted and simple(no self loop)\n\tedges will be tab (\\t) seperated (v1\tv2)\n";
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
	parse_args(argc, argv);//get input from command line
	graph_* g=new graph_(datafile);//read graph from data file
	
	//output parameters and graph statistics(edges and vertices)
	cout<<"Input file: "<<datafile<<"\nFraction of edges sampled: "<<sample_fraction<<"\n";
	cout<<"total edge: "<<g->getEdgeCount()<<"\n";
	cout<<"total vertices: "<<g->getVertexCount()<<"\n";
	cout.flush();
	
	time_tracker tt;//keep track of time of execution of graphlet count
	tt.start();
	cout<<"graphlet counts:\ng1:g2:g3:g4:g5:..........:g29:\n";
	g->signature_count_app_agg(sample_fraction);//Do graphlet count
	tt.stop();
	cout << "time: "<<tt.print()<<" second(s)"<<endl;
	delete g;
}
