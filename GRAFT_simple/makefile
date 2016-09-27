GRAFT.out : StringTokenizer.o GRAFT.cpp  graph.h  random.o time_tracker.h
	g++ -O3 -g -o GRAFT.out StringTokenizer.o random.o GRAFT.cpp 

StringTokenizer.o : StringTokenizer.cpp StringTokenizer.h
	g++ -O3 -c -g StringTokenizer.cpp StringTokenizer.h
	
random.o: random.h  random.cpp
	g++ -O3 -g -c random.h random.cpp

clean:
	rm -f *.o
	rm -f *.out
	rm -f *.gch
	rm -r *.dSYM
	rm -f *~
