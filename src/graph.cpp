#include "graph.h"
adjalist::adjalist(const adjalist &other)
{
     //copy construction function
	this->vernum=other.vernum;
	adjanodelist=new list<adjanode>[other.vernum];
	memcpy(this->adjanodelist,other.adjanodelist,(this->vernum)*sizeof(list<adjanode>));
}

bool operator>=(const adjanode& le,const adjanode& ri)
{
	return (le.dist>=ri.dist);
}
bool operator>(const adjanode& le,const adjanode& ri)
{
	return le.dist>ri.dist;
}
bool operator<=(const adjanode& le,const adjanode& ri)
{
	return le.dist<=ri.dist;
}
bool operator<(const adjanode& le,const adjanode& ri)
{
	return le.dist<ri.dist;
}
adjalist& adjalist::operator =(const adjalist &other)
{
	if (this==&other)
	{
		return (*this);	
	}
	if (adjanodelist!=NULL)
	{
		delete [] adjanodelist;
	}
	this->vernum=other.vernum;
	adjanodelist=new list<adjanode>[other.vernum];
	memcpy(this->adjanodelist,other.adjanodelist,(this->vernum)*sizeof(list<adjanode>));
	return (*this);
}
adjalist::~adjalist()
{
    if (this->adjanodelist!=NULL)
    {
		delete[] adjanodelist;
		adjanodelist=NULL;
    }
}

// bool ConvertHEMeshToAdjalist(HEMesh* hemesh,adjalist& adlist)
// {
//     if (hemesh==NULL)
//     {
// 		return false;
//     }
// 	if (hemesh->vnum==0)
// 	{
// 		return false;
// 	}
// 	if(adlist.vernum!=0)
// 	{
// 		return false;
// 	}
//     adlist.vernum=hemesh->vnum;
//     adlist.adjanodelist=new list<adjanode>[adlist.vernum];
//     for (int i=0;i<hemesh->vnum;i++)
//     {
// 		int outedgeindex=hemesh->verStartHEs[i];
// 		int tick=outedgeindex;
// 		bool boundary=false;
// 		do{
// 			int prev=hemesh->hedges[tick].prev;
// 			tick=hemesh->hedges[prev].opp;
// 			if(tick<0)
// 			{
// 				boundary=true;
// 				break;
// 			}
// 			//
// 			adjanode temp;
// 			temp.index=hemesh->hedges[tick].tv;
// 			int fvindex=hemesh->hedges[tick].fv;
// 			int tvindex=hemesh->hedges[tick].tv;
// 			temp.dist=hemesh->vertices[fvindex].Distance(hemesh->vertices[tvindex]);
// 			adlist.adjanodelist[i].push_back(temp);
// 		}while(tick!=outedgeindex);
// 		if(boundary)
// 		{
// 			tick=outedgeindex;
// 			do{
// 				int opp=hemesh->hedges[tick].opp;
// 				if(opp<0)
// 				{
// 					break;
// 				}
// 				tick=hemesh->hedges[opp].next;
// 				if(tick<0)
// 				{
// 					break;
// 				}
// 				adjanode temp;
// 				temp.index=hemesh->hedges[tick].tv;
// 				int fvindex=hemesh->hedges[tick].fv;
// 				int tvindex=hemesh->hedges[tick].tv;
// 				temp.dist=hemesh->vertices[fvindex].Distance(hemesh->vertices[tvindex]);
// 				adlist.adjanodelist[i].push_back(temp);
//    			}while(tick!=outedgeindex);
// 		}
//     }
// 	return true;
// }
bool ConvertDTriMeshToAdjalist(DTriMesh& mesh, adjalist& adlist)
{
	if (mesh.n_vertices()==0)
	{
		return false;
	}

    adlist.vernum=mesh.n_vertices();
    adlist.adjanodelist=new list<adjanode>[adlist.vernum];
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
	    DTriMesh::VertexHandle vi(i);
		OpenMesh::Vec3d pi = mesh.point(vi);
		for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi);vvi;vvi++)
		{
			OpenMesh::Vec3d pj = mesh.point(vvi);
			double _len = (pi-pj).length();
			adjanode temp;
			temp.index = vvi.handle().idx();
			temp.dist = _len;
			adlist.adjanodelist[i].push_back(temp);
		}
	}
	return true;
}



 GraPath::GraPath(const GraPath& other)
 {
    vernum=other.vernum;
	pathlist=new int[other.vernum];
	distlist=new double[other.vernum];
	memcpy(pathlist,other.pathlist,(other.vernum)*sizeof(int));
	memcpy(distlist,other.distlist,(other.vernum)*sizeof(double));
 };
GraPath& GraPath::operator=(const GraPath& other)
{
   //检查自赋值
   if (this==&other)
   {
	   return (*this);
   }
   //清空原有内存
   if (pathlist!=NULL)
   {
	   delete[] pathlist;
	   pathlist=NULL;
   }
   if (distlist!=NULL)
   {
	   delete[] distlist;
	   distlist=NULL;
   }
   //赋值
   this->vernum=other.vernum;
   pathlist=new int[other.vernum];
   memcpy(pathlist,other.pathlist,(other.vernum)*sizeof(int));
   //
   distlist=new double[other.vernum];
   memcpy(distlist,other.distlist,(other.vernum)*sizeof(double));
}
GraPath::~GraPath()
{
	if (pathlist!=NULL)
	{
		delete[] pathlist;
		pathlist=NULL;
	}
	if (distlist!=NULL)
	{
		delete[] distlist;
		distlist=NULL;
	}
}
bool operator<=(const DijNode& le,const DijNode& ri)
{
	return le.dist<=ri.dist;
}
bool operator<(const DijNode& le,const DijNode& ri)
{
	return le.dist<ri.dist;
}
bool operator>=(const DijNode& le,const DijNode& ri)
{
	return le.dist>=ri.dist;
}
bool operator>(const DijNode& le,const DijNode& ri)
{
	return le.dist>ri.dist;
}
bool Dijkstra(adjalist& meshadja,int source,GraPath& meshpath,int tar/* =NIL */)
{
	 if (source<0 || source>meshadja.vernum)
	 {
		 return false;
	 }
     double* dist=new double[meshadja.vernum];
	 int* previous=new int[meshadja.vernum];
	 int* vertimes=new int[meshadja.vernum];
	 //Initialization
	 for (int i=0;i<meshadja.vernum;i++)
	 {
		 dist[i]=DInf;
		 previous[i]=NIL;
		 vertimes[i]=0;
	 }
	 //
	 meshpath.vernum=meshadja.vernum;
	 meshpath.source=source;
	 meshpath.pathlist=new int[meshadja.vernum];
	 memset(meshpath.pathlist,NIL,meshadja.vernum*sizeof(int));
	 //
	 if (meshpath.distlist!=NULL)
	 {
		 delete[] meshpath.distlist;
	 }
	 //initialize
     meshpath.distlist=new double[meshadja.vernum];
	 for (int i=0;i<meshadja.vernum;i++)
	 {
		 meshpath.distlist[i]=DInf;
	 }
	 //from source to source
	 dist[source]=0;
//	 vertimes[source]=0;
	 Min_Heap<DijNode> Qset;
     DijNode tempnode(source,0,vertimes[source]);
	 Qset.Insert(tempnode);
	 while (!Qset.IsEmpty())
	 {
        DijNode u=Qset.Pop();
		if (u.times==vertimes[u.verindex])
		{
              if (u.dist==DInf)
              {
				  break;
              }
			  if (u.verindex==tar)
			  {
				  break;
			  }
			  for (list<adjanode>::iterator iter=meshadja.adjanodelist[u.verindex].begin();iter!=meshadja.adjanodelist[u.verindex].end();iter++)
              {
                   double alt=dist[u.verindex]+iter->dist;
				   if (alt<dist[iter->index])
				   {
					   dist[iter->index]=alt;
                       meshpath.pathlist[iter->index]=u.verindex;
                       vertimes[iter->index]++;
					   DijNode temp;
					   temp.dist=alt;
					   temp.verindex=iter->index;
					   temp.times=vertimes[iter->index];
					   Qset.Insert(temp);
				   }
              }
		}
	 } 
     //
	 memcpy(meshpath.distlist,dist,meshadja.vernum*sizeof(double));
	 //
// 	 for (int i=0;i<meshadja.vernum;i++)
// 	 {
// 	  	 meshpath.distlist[i]=dist[i];
// 	 }
	 // 
	 delete[] vertimes;
	 delete[] previous;
	 delete[] dist;
	 return true;
}

bool GraPath::GetGraPath(int tar, std::list<int> &grapath)
{
	int u=tar;
	if (u==NIL)
	{
		return false;
	}
	while (pathlist[u]!=NIL)
	{
		grapath.push_front(u);
		u=pathlist[u];
	}
	grapath.push_front(u);
	return true;
}

// bool Vertex2Edge(HEMesh& hemesh,list<int>& boundaryver,list<int>& boundaryedge)
// {
//      if (hemesh.vnum==0)
//      {
// 		 return false;
//      }
//      if (!boundaryedge.empty())
//      {
// 		 boundaryedge.clear();
//      }
// 	 if (boundaryver.empty())
// 	 {
// 		 return false;
// 	 }
// 	 list<int>::iterator iter=boundaryver.begin();
// 	 list<int>::iterator jiter=iter;
// 	 iter++;
// 	 while(iter!=boundaryver.end())
// 	 {
//          int heedge=GetHEEdge(hemesh,(*jiter),(*iter));
// 		 if (heedge!=(-1))
// 		 {
// 			 boundaryedge.push_back(heedge);
// 		 }
// 		 jiter=iter;
// 		 iter++;
// 	 }
// 	 return true;
// }