GRAFT_stratified.out : StringTokenizer.o GRAFT.cpp  graph.h  random.o triangle.h time_tracker.h
	g++ -O3 -g -o GRAFT_stratified.out -I /Users/mmrahman/Research/eclipse_cpp/boost_1_46_1 StringTokenizer.o random.o GRAFT.cpp 

StringTokenizer.o : StringTokenizer.cpp StringTokenizer.h
	g++ -O3 -c -g StringTokenizer.cpp StringTokenizer.h
	
random.o: random.h  random.cpp
	g++ -O3 -g -c -I /Users/mmrahman/Research/eclipse_cpp/boost_1_46_1 random.h random.cpp

clean:
	rm -f *.o
	rm -f *.out
	rm -f *.gch
	rm -r *.dSYM
	rm -f *~