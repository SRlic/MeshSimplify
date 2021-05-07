#include<fstream>
#include<sstream>
#include<string.h>
#include<iostream>
#include<pthread.h>
#include <iostream>
#include<pthread.h>
#include<unistd.h>
#include<time.h>
#include <windows.h>
#include"common.h"
using namespace std;

#define thread_num 16
#define thread_count 16
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
int Facenum, DelFacenum;
double Xmax,Xmin;
int findpart(double x,double *bound)
{
	for(int i=0;i<thread_num;i++)
	{
		if(x>=bound[i]&&x<=bound[i+1])
		{
			return i;
		}
	}

}
//按照x轴划分区域，标识边界边和点
void partition(E_Heap* eheap,V_Group* vGroup,init_message* message)
{
    cout<<"start partition"<<endl;
	double length=Xmax-Xmin;
	double avelen=length/thread_num;
	double bound[thread_num+1];
	for(int i=0;i<thread_num+1;i++)
	{
		bound[i]=Xmin+i*avelen;
	}
	for(int i=0;i<eheap->Edgelist.size();i++)
	{
	   // cout<<i<<endl;
	   // cout<<i<<endl;
		Edge e=eheap->Edgelist[i];
		if(findpart(vGroup->Vgroup[e.v1].pos.x,bound)!=findpart(vGroup->Vgroup[e.v2].pos.x,bound))
		{
            eheap->Edgelist[i].isbound=1;
            vGroup->Vgroup[e.v1].isbound=1;
            vGroup->Vgroup[e.v2].isbound=1;

		}
        else
        {
            int j=findpart(vGroup->Vgroup[e.v1].pos.x,bound);
            //cout<<j<<endl;
            message[j].idlist.push_back(e);
        }


	}
}
// E_Heap* eheap;
// V_Group* vGroup;
//void setRation(double ratio);
//void init(string path,E_Heap* eheap);
//void initall(E_Heap* eheap,V_Group* vGroup);
//void start(double ratio,E_Heap* eheap,V_Group* vGroup);
//void output(string path);
////计算顶点周围面误差
//Matrix calVertexDelta(int id，V_Group* vGroup);
//S_Vec3d calSimpPos(Edge e, Matrix m，V_Group* vGroup);
////边误差初始化
//void * P_Edge_init(void* args);
//void Edge_init(E_Heap* eheap,V_Group* vGroup);










Matrix calVertexDelta(int id,V_Group* vGroup)
{
	Matrix ans;
	Vertex *p = vGroup->Vgroup+id;
	for (set<int>::iterator it1 = p->Approachlist.begin(); it1 != p->Approachlist.end(); it1++) {
		for (set<int>::iterator it2 = p->Approachlist.begin(); it2 != p->Approachlist.end(); it2++) {
			if (*it1 < *it2 && NearTo(vGroup->Vgroup + (*it1), vGroup->Vgroup + (*it2)))
			{
				Vertex* v1=vGroup->Vgroup + (*it1);
				Vertex* v2 = vGroup->Vgroup + (*it2);
				Vec4 tempface = ComputeFace(p, v1, v2);
				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						ans.M[i][j] += tempface.M[i] * tempface.M[j];
					}
				}
			}
		}
	}
	return ans;
}

S_Vec3d calSimpPos(Edge e, Matrix m,V_Group *vGroup )
{
	m.M[3][0] = 0;
	m.M[3][1] = 0;
	m.M[3][2] = 0;
	m.M[3][3] = 1;

	Vec4 Y(0, 0, 0, 1);
	Vec4 ans = gass(m, Y);
	if (ans.M[3] > EPS)
	{
		return S_Vec3d(ans.M[0], ans.M[1], ans.M[2]);
	}
	else
	{
		return S_Vec3d((vGroup->Vgroup[e.v1].pos.x + vGroup->Vgroup[e.v2].pos.x) / 2,
			(vGroup->Vgroup[e.v1].pos.y + vGroup->Vgroup[e.v2].pos.y) / 2,
			(vGroup->Vgroup[e.v1].pos.z + vGroup->Vgroup[e.v2].pos.z) / 2
		);
	}
}


void readfile(string path,E_Heap* eheap,V_Group* vGroup)
{
    cout<<"start read file"<<endl;
	ifstream file(path);
	string line;
	while (getline(file, line))
	{
		if (line.substr(0, 2) == "vt")
		{

		}
		else if (line.substr(0, 2) == "vn")
		{

		}
		else if (line.substr(0, 1) == "v")
		{
			//更新点集
			double x, y, z;
			istringstream s(line.substr(2));
			s >> x; s >> y; s >> z;
			if(x>Xmax)
			{
				Xmax=x;
			}
			if(x<Xmin)
			{
				Xmin=x;
			}
			vGroup->add_V(Vertex(x, y, z));

		}
		else if (line.substr(0, 1) == "f")
		{
			//更新点的邻接表以及边集
			int a, b, c;
			istringstream vtns(line.substr(2));
			vtns >> a; vtns >> b; vtns >> c;
			Ver_insert(vGroup->Vgroup + a, b);
			Ver_insert(vGroup->Vgroup + a, c);
			Ver_insert(vGroup->Vgroup + b, a);
			Ver_insert(vGroup->Vgroup + b, c);
			Ver_insert(vGroup->Vgroup + c, a);
			Ver_insert(vGroup->Vgroup + c, b);


			Edge e1(a, b);
			if(eheap->mapEdgeToID.count(make_pair(e1.v1, e1.v2))<=0)
			{
				eheap->addEdge(e1);
			}
			Edge e2(b, c);
			if (eheap->mapEdgeToID.count(make_pair(e2.v1, e2.v2)) <= 0)
			{
				eheap->addEdge(e2);
			}
			Edge e3(a, c);
			if (eheap->mapEdgeToID.count(make_pair(e3.v1, e3.v2)) <= 0)
			{
				eheap->addEdge(e3);
			}
			Facenum++;

		}
		else if (line.substr(0, 1) == "#")
		{

		}
		else
		{

		}

	}
	file.close();
}

Edge getMinD(priority_queue<Edge, std::vector<Edge>, E_cmp> EdgeQueue,E_Heap* eheap)
{
	if (EdgeQueue.size() <= 0) {
		return Edge(0, 0);
	}
	while (eheap->isDeleted[EdgeQueue.top().id]) {
		EdgeQueue.pop();
	}
	Edge e = EdgeQueue.top();
	EdgeQueue.pop();
	return e;
	return Edge();

}



void * P_meshsimple(void* args)
{
//    cout<<"ssss'"<<endl;
	init_message *message;
	message=(init_message *) args;
//	int start=message->start;
	//int end=message->end;
	int av_Del=message->av_Del;
	E_Heap* eheap=message->eheap;
	V_Group* vGroup=message->vGroup;
	//初始化计算误差，并将不是边界边排序
    //pthread_mutex_lock(&mutex);
	for (int i=0;i<message->idlist.size();i++)
	{
		Edge e1 = message->idlist[i];
		Matrix mat = calVertexDelta(e1.v1,vGroup) + calVertexDelta(e1.v2,vGroup);
		e1.simp_v = calSimpPos(e1,mat,vGroup);
		Vec4 X(e1.simp_v.x, e1.simp_v.y, e1.simp_v.z, 1.0);
		if (vGroup->getCommonVertexNum(e1.v1, e1.v2) != 2) {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			e1.deItav = 0;
		//	cout<<"error quit"<<endl;
		//	pthread_exit(NULL);
		}
		double pri = 0;
		for (int i = 0; i < 4; i++) {
			double p = 0;
			for (int j = 0; j < 4; j++)
				p += X.M[j] * mat.M[i][j];
			pri += p * X.M[i];
		}
		e1.deItav = pri;
		//message->result.push_back(e1);
		if(e1.isbound==0)
		{
			message->EdgeQueue.push(e1);
		}
		//eheap.EdgeQueue.push(e1);
	}
	//pthread_mutex_unlock(&mutex);

	//开始删边

	// cout<<av_Del<<endl;
	// av_Del=0; // 跳过迭代过程
	 for(int i=0;i<av_Del;i+=2)
	 {
       // cout<<i<<endl;

		Edge e = getMinD(message->EdgeQueue,eheap);
		Vertex* v1 = &(vGroup->Vgroup[e.v1]);
		Vertex* v2 = &(vGroup->Vgroup[e.v2]);
		//新顶点
		Vertex  tempv(e.simp_v.x, e.simp_v.y,e.simp_v.z);
		pthread_mutex_lock(&mutex);
		int v0_id = vGroup->add_V(tempv);
		pthread_mutex_unlock(&mutex);
		Vertex* v0 = vGroup->Vgroup + v0_id;

		set<int> approachlist0;
		approachlist0.clear();
		//删除旧边
		pthread_mutex_lock(&mutex);
		eheap->delEdge(e);
		//更新新点/旧点临近点
		for (set<int>::iterator it = v1->Approachlist.begin(); it != v1->Approachlist.end(); it++) {
			if ((*it) != v2->id) {
				//删除旧边
				eheap->delEdge(Edge((*it), v1->id));
				vGroup->Vgroup[(*it)].Approachlist.erase(v1->id);
				approachlist0.insert((*it));
			}

		}
		for (set<int>::iterator it = v2->Approachlist.begin(); it != v2->Approachlist.end(); it++) {
			if ((*it) != v1->id) {
				//删除旧边
				eheap->delEdge(Edge((*it), v2->id));
				vGroup->Vgroup[(*it)].Approachlist.erase(v2->id);
				approachlist0.insert((*it));
			}

		}
		for (set<int>::iterator it = approachlist0.begin(); it != approachlist0.end(); it++) {
			vGroup->Vgroup[(*it)].Approachlist.insert(v0_id);
			vGroup->Vgroup[v0_id].Approachlist.insert(*it);
		}
		//删除顶点
		pthread_mutex_unlock(&mutex);
		vGroup->isDeleted[v1->id] = true;
		vGroup->isDeleted[v2->id] = true;

		//更新边堆
		for (set<int>::iterator it = approachlist0.begin(); it != approachlist0.end(); it++) {
			Edge e1((*it), v0_id);
				pthread_mutex_lock(&mutex);
			eheap->addEdge(e1);
					pthread_mutex_unlock(&mutex);
            Matrix mat = calVertexDelta(e1.v1,vGroup) + calVertexDelta(e1.v2,vGroup);
            e1.simp_v = calSimpPos(e1,mat,vGroup);
			Vec4 X(e1.simp_v.x, e1.simp_v.y, e1.simp_v.z, 1.0);
			if (vGroup->getCommonVertexNum(e1.v1, e1.v2) != 2) {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				e1.deItav = 0;
				//return;
			}
			double pri = 0;
			for (int i = 0; i < 4; i++) {
				double p = 0;
				for (int j = 0; j < 4; j++)
					p += X.M[j] * mat.M[i][j];
				pri += p * X.M[i];
			}
			e1.deItav = pri;
			message->EdgeQueue.push(e1);
		}
	 }

}

void P_start(double ratio,E_Heap* eheap,V_Group* vGroup)
{
    cout<<"start_init"<<endl;
	long thread;
	pthread_t* thread_handles;
	//均分删边数
	DelFacenum = Facenum * (1 - ratio); //计算删除面的数量
	cout<<DelFacenum<<endl;
	int av_Del=DelFacenum/thread_count;
	//均分边
	thread_handles=(pthread_t*)  malloc(thread_count*sizeof(pthread_t));
    int edge_count=eheap->Edgelist.size();
	int p_num=edge_count/thread_count+1;
	//消息初始化V_Group()
	init_message message[thread_count];
	for(int i=0;i<thread_count;i++)
	{
		//message[i].start=i*p_num;
		//message[i].end=(i+1)*p_num;
		//message[i].av_Del=av_Del;
//		message[i].eheap=new E_Heap();
//		message[i].eheap->Edgelist=eheap->Edgelist;
//		message[i].eheap->mapEdgeToID=eheap->mapEdgeToID;
//		message[i].vGroup=new V_Group();
//		for(int j=0;j<MAX_SIZE;j++)
//        {
//            message[i].vGroup->Vgroup[j]=vGroup->Vgroup[j];
//        }
        message[i].eheap=eheap;
        message[i].vGroup=vGroup;

	}
	 partition(eheap,vGroup,message);
    for(int i=0;i<thread_count;i++)
	{
	    message[i].av_Del=(message[i].idlist.size()*1.0/eheap->Edgelist.size()*1.0)*DelFacenum;
	}
	//开始删边
	for(thread=0;thread<thread_count;thread++)
	{

		pthread_create(&thread_handles[thread],NULL,P_meshsimple,(void*) &message[thread]);

	}
	for(thread=0;thread<thread_count;thread++)
	{
		pthread_join(thread_handles[thread],NULL);

	}

	free(thread_handles);
}
void * P_meshsimple_random1(void* args)
{
	init_message *message;
	message=(init_message *) args;
	int start=message->start;
	int end=message->end;
	int av_Del=message->av_Del;
	E_Heap* eheap=message->eheap;
	V_Group* vGroup=message->vGroup;
	//初始化计算误差，并将不是边界边排序
    //pthread_mutex_lock(&mutex);
	for (int i=start;i<end&&i<eheap->Edgelist.size();i++)
	{
		Edge e1 = eheap->Edgelist[i];
		Matrix mat = calVertexDelta(e1.v1,vGroup) + calVertexDelta(e1.v2,vGroup);
		e1.simp_v = calSimpPos(e1,mat,vGroup);
		Vec4 X(e1.simp_v.x, e1.simp_v.y, e1.simp_v.z, 1.0);
		if (vGroup->getCommonVertexNum(e1.v1, e1.v2) != 2) {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			e1.deItav = 0;
		//	cout<<"error quit"<<endl;
		//	pthread_exit(NULL);
		}
		double pri = 0;
		for (int i = 0; i < 4; i++) {
			double p = 0;
			for (int j = 0; j < 4; j++)
				p += X.M[j] * mat.M[i][j];
			pri += p * X.M[i];
		}
		e1.deItav = pri;
		//message->result.push_back(e1);
		if(e1.isbound==0)
		{
			message->EdgeQueue.push(e1);
		}
		//eheap.EdgeQueue.push(e1);
	}
	//pthread_mutex_unlock(&mutex);

	//开始删边

	// cout<<av_Del<<endl;
	// av_Del=0; // 跳过迭代过程
	 for(int i=0;i<av_Del;i+=2)
	 {
       // cout<<i<<endl;

		Edge e = getMinD(message->EdgeQueue,eheap);
		Vertex* v1 = &(vGroup->Vgroup[e.v1]);
		Vertex* v2 = &(vGroup->Vgroup[e.v2]);
		//新顶点
		Vertex  tempv(e.simp_v.x, e.simp_v.y,e.simp_v.z);
		pthread_mutex_lock(&mutex);
		int v0_id = vGroup->add_V(tempv);
		pthread_mutex_unlock(&mutex);
		Vertex* v0 = vGroup->Vgroup + v0_id;

		set<int> approachlist0;
		approachlist0.clear();
		//删除旧边
		pthread_mutex_lock(&mutex);
		eheap->delEdge(e);
		//更新新点/旧点临近点
		for (set<int>::iterator it = v1->Approachlist.begin(); it != v1->Approachlist.end(); it++) {
			if ((*it) != v2->id) {
				//删除旧边
				eheap->delEdge(Edge((*it), v1->id));
				vGroup->Vgroup[(*it)].Approachlist.erase(v1->id);
				approachlist0.insert((*it));
			}

		}
		for (set<int>::iterator it = v2->Approachlist.begin(); it != v2->Approachlist.end(); it++) {
			if ((*it) != v1->id) {
				//删除旧边
				eheap->delEdge(Edge((*it), v2->id));
				vGroup->Vgroup[(*it)].Approachlist.erase(v2->id);
				approachlist0.insert((*it));
			}

		}
		for (set<int>::iterator it = approachlist0.begin(); it != approachlist0.end(); it++) {
			vGroup->Vgroup[(*it)].Approachlist.insert(v0_id);
			vGroup->Vgroup[v0_id].Approachlist.insert(*it);
		}
		//删除顶点
		pthread_mutex_unlock(&mutex);
		vGroup->isDeleted[v1->id] = true;
		vGroup->isDeleted[v2->id] = true;

		//更新边堆
		for (set<int>::iterator it = approachlist0.begin(); it != approachlist0.end(); it++) {
			Edge e1((*it), v0_id);
				pthread_mutex_lock(&mutex);
			eheap->addEdge(e1);
					pthread_mutex_unlock(&mutex);
            Matrix mat = calVertexDelta(e1.v1,vGroup) + calVertexDelta(e1.v2,vGroup);
            e1.simp_v = calSimpPos(e1,mat,vGroup);
			Vec4 X(e1.simp_v.x, e1.simp_v.y, e1.simp_v.z, 1.0);
			if (vGroup->getCommonVertexNum(e1.v1, e1.v2) != 2) {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				e1.deItav = 0;
				//return;
			}
			double pri = 0;
			for (int i = 0; i < 4; i++) {
				double p = 0;
				for (int j = 0; j < 4; j++)
					p += X.M[j] * mat.M[i][j];
				pri += p * X.M[i];
			}
			e1.deItav = pri;
			message->EdgeQueue.push(e1);
		}
	 }

}

void P_start_random1(double ratio,E_Heap* eheap,V_Group* vGroup)
{
    cout<<"start_init"<<endl;
	long thread;
	pthread_t* thread_handles;
	//均分删边数
	DelFacenum = Facenum * (1 - ratio); //计算删除面的数量
	cout<<DelFacenum<<endl;
	int av_Del=DelFacenum/thread_count;
	//均分边
	thread_handles=(pthread_t*)  malloc(thread_count*sizeof(pthread_t));
    int edge_count=eheap->Edgelist.size();
	int p_num=edge_count/thread_count+1;
	//消息初始化V_Group()
	init_message message[thread_count];
	for(int i=0;i<thread_count;i++)
	{
		message[i].start=i*p_num;
		message[i].end=(i+1)*p_num;
		message[i].av_Del=av_Del;
//		message[i].eheap=new E_Heap();
//		message[i].eheap->Edgelist=eheap->Edgelist;
//		message[i].eheap->mapEdgeToID=eheap->mapEdgeToID;
//		message[i].vGroup=new V_Group();
//		for(int j=0;j<MAX_SIZE;j++)
//        {
//            message[i].vGroup->Vgroup[j]=vGroup->Vgroup[j];
//        }
        message[i].eheap=eheap;
        message[i].vGroup=vGroup;

	}
	//partition(eheap,vGroup,message);
    // for(int i=0;i<thread_count;i++)
	// {
	//     message[i].av_Del=(message[i].idlist.size()*1.0/eheap->Edgelist.size()*1.0)*DelFacenum;
	// }
	//开始删边
	for(thread=0;thread<thread_count;thread++)
	{

		pthread_create(&thread_handles[thread],NULL,P_meshsimple_random1,(void*) &message[thread]);

	}
	for(thread=0;thread<thread_count;thread++)
	{
		pthread_join(thread_handles[thread],NULL);

	}

	free(thread_handles);
}

void * P_meshsimple_random2(void* args)
{
	init_message *message;
	message=(init_message *) args;
	// int start=message->start;
	// int end=message->end;
	int av_Del=message->av_Del;
	int split=message->split;
	E_Heap* eheap=message->eheap;
	V_Group* vGroup=message->vGroup;
	//初始化计算误差，并将不是边界边排序
    //pthread_mutex_lock(&mutex);
	for (int i=0;i<eheap->Edgelist.size();i+=split)
	{
		Edge e1 = eheap->Edgelist[i];
		Matrix mat = calVertexDelta(e1.v1,vGroup) + calVertexDelta(e1.v2,vGroup);
		e1.simp_v = calSimpPos(e1,mat,vGroup);
		Vec4 X(e1.simp_v.x, e1.simp_v.y, e1.simp_v.z, 1.0);
		if (vGroup->getCommonVertexNum(e1.v1, e1.v2) != 2) {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			e1.deItav = 0;
		//	cout<<"error quit"<<endl;
		//	pthread_exit(NULL);
		}
		double pri = 0;
		for (int i = 0; i < 4; i++) {
			double p = 0;
			for (int j = 0; j < 4; j++)
				p += X.M[j] * mat.M[i][j];
			pri += p * X.M[i];
		}
		e1.deItav = pri;
		//message->result.push_back(e1);
		if(e1.isbound==0)
		{
			message->EdgeQueue.push(e1);
		}
		//eheap.EdgeQueue.push(e1);
	}
	//pthread_mutex_unlock(&mutex);

	//开始删边

	// cout<<av_Del<<endl;
	// av_Del=0; // 跳过迭代过程
	 for(int i=0;i<av_Del;i+=2)
	 {
       // cout<<i<<endl;

		Edge e = getMinD(message->EdgeQueue,eheap);
		Vertex* v1 = &(vGroup->Vgroup[e.v1]);
		Vertex* v2 = &(vGroup->Vgroup[e.v2]);
		//新顶点
		Vertex  tempv(e.simp_v.x, e.simp_v.y,e.simp_v.z);
		pthread_mutex_lock(&mutex);
		int v0_id = vGroup->add_V(tempv);
		pthread_mutex_unlock(&mutex);
		Vertex* v0 = vGroup->Vgroup + v0_id;

		set<int> approachlist0;
		approachlist0.clear();
		//删除旧边
		pthread_mutex_lock(&mutex);
		eheap->delEdge(e);
		//更新新点/旧点临近点
		for (set<int>::iterator it = v1->Approachlist.begin(); it != v1->Approachlist.end(); it++) {
			if ((*it) != v2->id) {
				//删除旧边
				eheap->delEdge(Edge((*it), v1->id));
				vGroup->Vgroup[(*it)].Approachlist.erase(v1->id);
				approachlist0.insert((*it));
			}

		}
		for (set<int>::iterator it = v2->Approachlist.begin(); it != v2->Approachlist.end(); it++) {
			if ((*it) != v1->id) {
				//删除旧边
				eheap->delEdge(Edge((*it), v2->id));
				vGroup->Vgroup[(*it)].Approachlist.erase(v2->id);
				approachlist0.insert((*it));
			}

		}
		for (set<int>::iterator it = approachlist0.begin(); it != approachlist0.end(); it++) {
			vGroup->Vgroup[(*it)].Approachlist.insert(v0_id);
			vGroup->Vgroup[v0_id].Approachlist.insert(*it);
		}
		//删除顶点
		pthread_mutex_unlock(&mutex);
		vGroup->isDeleted[v1->id] = true;
		vGroup->isDeleted[v2->id] = true;

		//更新边堆
		for (set<int>::iterator it = approachlist0.begin(); it != approachlist0.end(); it++) {
			Edge e1((*it), v0_id);
				pthread_mutex_lock(&mutex);
			eheap->addEdge(e1);
					pthread_mutex_unlock(&mutex);
            Matrix mat = calVertexDelta(e1.v1,vGroup) + calVertexDelta(e1.v2,vGroup);
            e1.simp_v = calSimpPos(e1,mat,vGroup);
			Vec4 X(e1.simp_v.x, e1.simp_v.y, e1.simp_v.z, 1.0);
			if (vGroup->getCommonVertexNum(e1.v1, e1.v2) != 2) {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				e1.deItav = 0;
				//return;
			}
			double pri = 0;
			for (int i = 0; i < 4; i++) {
				double p = 0;
				for (int j = 0; j < 4; j++)
					p += X.M[j] * mat.M[i][j];
				pri += p * X.M[i];
			}
			e1.deItav = pri;
			message->EdgeQueue.push(e1);
		}
	 }

}
void P_start_random2(double ratio,E_Heap* eheap,V_Group* vGroup)
{
    cout<<"start_init"<<endl;
	long thread;
	pthread_t* thread_handles;
	//均分删边数
	DelFacenum = Facenum * (1 - ratio); //计算删除面的数量
	cout<<DelFacenum<<endl;
	int av_Del=DelFacenum/thread_count;
	//均分边
	thread_handles=(pthread_t*)  malloc(thread_count*sizeof(pthread_t));
    int edge_count=eheap->Edgelist.size();
	int p_num=edge_count/thread_count+1;
	//消息初始化V_Group()
	init_message message[thread_count];
	for(int i=0;i<thread_count;i++)
	{
		message[i].start=i*p_num;
		message[i].end=(i+1)*p_num;
		message[i].av_Del=av_Del;
		message[i].split=thread_count;
//		message[i].eheap=new E_Heap();
//		message[i].eheap->Edgelist=eheap->Edgelist;
//		message[i].eheap->mapEdgeToID=eheap->mapEdgeToID;
//		message[i].vGroup=new V_Group();
//		for(int j=0;j<MAX_SIZE;j++)
//        {
//            message[i].vGroup->Vgroup[j]=vGroup->Vgroup[j];
//        }
        message[i].eheap=eheap;
        message[i].vGroup=vGroup;

	}
	//partition(eheap,vGroup,message);
    // for(int i=0;i<thread_count;i++)
	// {
	//     message[i].av_Del=(message[i].idlist.size()*1.0/eheap->Edgelist.size()*1.0)*DelFacenum;
	// }
	//开始删边
	for(thread=0;thread<thread_count;thread++)
	{

		pthread_create(&thread_handles[thread],NULL,P_meshsimple_random2,(void*) &message[thread]);

	}
	for(thread=0;thread<thread_count;thread++)
	{
		pthread_join(thread_handles[thread],NULL);

	}

	free(thread_handles);
}


void start(double ratio,E_Heap* eheap,V_Group* vGroup)
{
	//Edge_init();  //初始化边队列
	DelFacenum = Facenum * (1 - ratio); //计算删除面的数量
    cout<<"start simplfy"<<endl;
    cout<<DelFacenum<<endl;
   // cout<<Facenum<<endl;
	for (int i = 0; i < DelFacenum; i += 2)
	{
		Edge e = eheap->getMinD();
		Vertex* v1 = &(vGroup->Vgroup[e.v1]);
		Vertex* v2 = &(vGroup->Vgroup[e.v2]);
		//新顶点
		Vertex  tempv(e.simp_v.x, e.simp_v.y,e.simp_v.z);
		int v0_id = vGroup->add_V(tempv);
		Vertex* v0 = vGroup->Vgroup + v0_id;

		set<int> approachlist0;
		approachlist0.clear();
		//删除旧边
		eheap->delEdge(e);
		//更新新点/旧点临近点
		for (set<int>::iterator it = v1->Approachlist.begin(); it != v1->Approachlist.end(); it++) {
			if ((*it) != v2->id) {
				//删除旧边
				eheap->delEdge(Edge((*it), v1->id));
				vGroup->Vgroup[(*it)].Approachlist.erase(v1->id);
				approachlist0.insert((*it));
			}

		}

		for (set<int>::iterator it = v2->Approachlist.begin(); it != v2->Approachlist.end(); it++) {
			if ((*it) != v1->id) {
				//删除旧边
				eheap->delEdge(Edge((*it), v2->id));
				vGroup->Vgroup[(*it)].Approachlist.erase(v2->id);
				approachlist0.insert((*it));
			}

		}
		for (set<int>::iterator it = approachlist0.begin(); it != approachlist0.end(); it++) {
			vGroup->Vgroup[(*it)].Approachlist.insert(v0_id);
			vGroup->Vgroup[v0_id].Approachlist.insert(*it);
		}
		//删除顶点

		vGroup->isDeleted[v1->id] = true;
		vGroup->isDeleted[v2->id] = true;

		//更新边堆
		for (set<int>::iterator it = approachlist0.begin(); it != approachlist0.end(); it++) {
			Edge e1((*it), v0_id);
			eheap->addEdge(e1);
            Matrix mat = calVertexDelta(e1.v1,vGroup) + calVertexDelta(e1.v2,vGroup);
            e1.simp_v = calSimpPos(e1,mat,vGroup);
			Vec4 X(e1.simp_v.x, e1.simp_v.y, e1.simp_v.z, 1.0);
			if (vGroup->getCommonVertexNum(e1.v1, e1.v2) != 2) {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				e1.deItav = 0;
				//return;
			}
			double pri = 0;
			for (int i = 0; i < 4; i++) {
				double p = 0;
				for (int j = 0; j < 4; j++)
					p += X.M[j] * mat.M[i][j];
				pri += p * X.M[i];
			}
			e1.deItav = pri;
			eheap->EdgeQueue.push(e1);
		}
	}

}

void output(string path,E_Heap* eheap,V_Group* vGroup)
{
	map<int, int> oldTonew;
	ofstream file(path);
	int cnt = 0;
	int cntv = 0, cntf = 0;
	for (int i = 1; i <= vGroup->curnum; i++) {//输出所有点
		if (vGroup->isDeleted[i])//如果第i个点已经删掉了，就略去
			continue;
		Vertex* v = &vGroup->Vgroup[i];
		cnt++;
		oldTonew.insert(make_pair(v->id, cnt));
		file << "v" << " " << v->pos.x << " " << v->pos.y << " " << v->pos.z << endl;
	}
	for (int i = 1; i <= vGroup->curnum; i++) {//输出所有面
		if (vGroup->isDeleted[i])//如果第i个点已经删掉了，就略去
			continue;
		Vertex* v = &(vGroup->Vgroup[i]);//对于第i个点
		for (set<int>::iterator it1 = v->Approachlist.begin(); it1 != v->Approachlist.end(); it1++) {
			if (i >= (*it1))
				continue;
			for (set<int>::iterator it2 = v->Approachlist.begin(); it2 != v->Approachlist.end(); it2++) {
				if ((*it1) < (*it2) && NearTo(vGroup->Vgroup+(*it1), vGroup->Vgroup + (*it2)) )
				{
					file << "f" << " " << oldTonew[v->id] << " " << oldTonew[vGroup->Vgroup[(*it1)].id] << " " << oldTonew[vGroup->Vgroup[(*it2)].id] << endl;
					//printf("f %d %d %d\n", v->id, vGroup->Vgroup[(*it1)].id, vGroup->Vgroup[(*it2)].id);
					cntf++;
				}
			}
		}

	}

}

void Edge_init(E_Heap* eheap,V_Group* vGroup)
{
    cout<<"start init"<<endl;
	vector<Edge>::iterator it = eheap->Edgelist.begin();
	for (it = eheap->Edgelist.begin(); it != eheap->Edgelist.end(); it++)
	{
		Edge e1 = *it;
		//边各点移动误差计算以及坍塌后点确定
		Matrix mat = calVertexDelta(e1.v1,vGroup) + calVertexDelta(e1.v2,vGroup);
		e1.simp_v = calSimpPos(e1,mat,vGroup);
		Vec4 X(e1.simp_v.x, e1.simp_v.y, e1.simp_v.z, 1.0);
		if (vGroup->getCommonVertexNum(e1.v1, e1.v2) != 2) {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			e1.deItav = 0;
			//cout<<"error,quit"<<endl;
			//return;
		}
		double pri = 0;
		for (int i = 0; i < 4; i++) {
			double p = 0;
			for (int j = 0; j < 4; j++)
				p += X.M[j] * mat.M[i][j];
			pri += p * X.M[i];
		}
		//坍塌后代价
		e1.deItav = pri;
		//堆排序插入
		eheap->EdgeQueue.push(e1);
	}
}

void simple(double ratio,string s1,string s2,string s3)
{
    long long head, tail, freq;
	Xmax=-100000;
	Xmin=100000;
    Facenum=0;
    DelFacenum=0;
    E_Heap* eheap=new E_Heap();
    V_Group* vGroup=new V_Group();
    readfile(s1,eheap,vGroup);
//	partition(eheap,vGroup);
//	P_start(ratio,eheap,vGroup);
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);	// similar to CLOCKS_PER_SEC
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
	P_start(ratio,eheap,vGroup);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "simple1: " << (tail - head) * 1000.0 / freq << "ms" << endl;


    output(s3,eheap,vGroup);
// 	Facenum=0;
//    DelFacenum=0;
//    eheap=new E_Heap();
//    vGroup=new V_Group();
//    readfile(s1,eheap,vGroup);
//    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);	// similar to CLOCKS_PER_SEC
//    QueryPerformanceCounter((LARGE_INTEGER*)&head);
//    P_start_random1(ratio,eheap,vGroup);
//    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//    cout << "simple_random1: " << (tail - head) * 1000.0 / freq << "ms" << endl;


 //   output(s3,eheap,vGroup);
// 	Facenum=0;
//    DelFacenum=0;
//    eheap=new E_Heap();
//    vGroup=new V_Group();
//    readfile(s1,eheap,vGroup);
//    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);	// similar to CLOCKS_PER_SEC
//    QueryPerformanceCounter((LARGE_INTEGER*)&head);
//    P_start_random2(ratio,eheap,vGroup);
//    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//    cout << "simple_random2: " << (tail - head) * 1000.0 / freq << "ms" << endl;

    Facenum=0;
    DelFacenum=0;
    eheap=new E_Heap();
    vGroup=new V_Group();
    readfile(s1,eheap,vGroup);
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);	// similar to CLOCKS_PER_SEC
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    Edge_init(eheap,vGroup);
    start(ratio,eheap,vGroup);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "simple2: " << (tail - head) * 1000.0 / freq << "ms" << endl;
//    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);	// similar to CLOCKS_PER_SEC
//    QueryPerformanceCounter((LARGE_INTEGER*)&head);
//    start(ratio,eheap,vGroup);
//    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//    cout << "simple2: " << (tail - head) * 1000.0 / freq << "ms" << endl;
    output(s2,eheap,vGroup);
    //printf("myerror\n");



}

int main()
{
    simple(0.5,"D:\\9_19_test\\kitten.obj","D:\\9_19_test\\kitten50.obj","D:\\9_19_test\\kitten50_P.obj");
    printf("myerror\n");
}
