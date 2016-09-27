/*
 * triangle.h
 *
 *  Created on: Nov 23, 2011
 *      Author: mmrahman
 */

#ifndef TRIANGLE_H_
#define TRIANGLE_H_

using namespace std;

class triangle
{
private:
	vector<int> vertices;
	int sc;
public:
	triangle(int X,int Y,int Z)
	{
		sc=0;
		if(X==Y or Y==Z or Z==X){
			cout<<X<<","<<Y<<","<<Z<<"\n";
			cout<<"triangle parameter is wrong!!!\n";
			exit(1);
		}
		vertices.push_back(X);
		vertices.push_back(Y);
		vertices.push_back(Z);
		sort(vertices.begin(),vertices.end());
	}

	~triangle(){
		vertices.clear();
		//cout<<"triangle destructor\n";
	}

	void increase_sc(){
		sc++;
	}

	int get_sc(){
		return sc;
	}
	void getVertices(vector<int> & t)
	{
		t.push_back(vertices[0]);
		t.push_back(vertices[1]);
		t.push_back(vertices[2]);
	}

	void print_triangle(){
		cout<<"tri: ("<< vertices[0] << "," << vertices[1] << "," << vertices[2] << "): "<< sc <<"\n";
	}
};

#endif /* TRIANGLE_H_ */
