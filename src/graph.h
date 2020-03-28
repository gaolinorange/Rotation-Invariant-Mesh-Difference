#pragma once
//#include "../geometry/mesh.h"
#include <list>
#include "Heap.h"
#define NIL -1
#define Inf 2147483647
#define DInf 1e+16
#include "Align.h"
//using namespace Genie::Geometry;

using namespace std;
// graph 
class adjanode{
public:
	int index;
	double dist;
	adjanode():index(0),dist(0) {}
	friend bool operator>=(const adjanode& le,const adjanode& ri);
	friend bool operator>(const adjanode& le,const adjanode& ri);
	friend bool operator<=(const adjanode& le,const adjanode& ri);
	friend bool operator<(const adjanode& le,const adjanode& ri);
};

class adjalist{
public:
	int vernum;
    list<adjanode>* adjanodelist;
	adjalist():vernum(0),adjanodelist(NULL) {}
	adjalist(const adjalist& other);
	adjalist& operator=(const adjalist& other);    
	~adjalist();
};

class GraPath
{
public:
	int vernum;
	int* pathlist;
	double* distlist;
	int source;
	GraPath():vernum(0),pathlist(NULL),source(0),distlist(NULL) {}
    bool GetGraPath(int tar,list<int>& grapath);
	GraPath(int _vernum)
	{
		this->vernum=_vernum;
		pathlist=new int[_vernum];
		memset(pathlist,NIL,_vernum*sizeof(int));
		distlist=new double[_vernum];
		for (int i=0;i<_vernum;i++)
		{
			distlist[i]=DInf;
		}
	}
	GraPath(const GraPath& other);
	GraPath& operator=(const GraPath& other);
	~GraPath();
};

class DijNode
{
public:
	int verindex;
	double dist;
	int times;
	DijNode():verindex(NIL),dist(0),times(0) {}
	DijNode(int _verindex,double _dist, int _times) {verindex=_verindex;dist=_dist;times=_times;}
	friend bool operator<=(const DijNode& le,const DijNode& ri);
	friend bool operator<(const DijNode& le,const DijNode& ri);
	friend bool operator>=(const DijNode& le,const DijNode& ri);
	friend bool operator>(const DijNode& le,const DijNode& ri);
};

//bool ConvertHEMeshToAdjalist(HEMesh* hemesh,adjalist& adlist);

bool ConvertDTriMeshToAdjalist(DTriMesh& mesh, adjalist& adlist);

bool Dijkstra(adjalist& meshadja,int source,GraPath& meshpath,int tar=NIL);

//bool Vertex2Edge(HEMesh& hemesh,list<int>& boundaryver,list<int>& boundaryedge);

