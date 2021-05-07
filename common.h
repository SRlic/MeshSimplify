#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

/*
*     共用数据结构，用于lod和mesh simplify
*     三维向量 S_Vec3d /边 Edge /点 Vertex /四维矩阵 Matrix
*     点group  V_Group  /   边堆  E_Heap
*     Date: 2020-9-18
*
*/
#include<iostream>
#include<vector>
#include <set>
#include <map>
#include <queue>
#include<math.h>
using namespace std;



#define MAX_SIZE 1000*1024
const double EPS= 1e-8;
//三维向量
typedef struct S_Vec3d
{
	S_Vec3d()
	{
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}
	S_Vec3d(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	double x;
	double y;
	double z;

}S_Vec3d;

//边
typedef struct Edge
{
	//注意顶点编号v1<v2
	Edge()
	{
		this->id = -1;
		this->v1 = -1;
		this->v2 = -1;
		this->deItav = 1000000;

	}
	Edge(int a, int b)
	{
		this->id = -1;
		if (a < b)
		{
			this->v1 = a;
			this->v2 = b;
		}
		else
		{
			this->v1 = b;
			this->v2 = a;
		}
		this->deItav = 1000000; // 赋一个无穷大值
	}
	Edge(int a, int b, int _id)
	{
		this->id = _id;
		if (a < b)
		{
			this->v1 = a;
			this->v2 = b;
		}
		else
		{
			this->v1 = b;
			this->v2 = a;
		}
		this->deItav = 1000000; // 赋一个无穷大值
	}
	int id;//边标号
	int v1, v2; //顶点标号
	S_Vec3d simp_v; //边坍塌后点
	double deItav; //边坍塌代价
	int isbound;

}Edge;
//点
class Vertex
{
public:
	int id;
	S_Vec3d pos; //三维坐标
	set<int> Approachlist; //邻接点表
	int isbound;
	Vertex()
	{
		this->id = -1;
		this->pos = S_Vec3d();
		//this->Approachlist.clear();
		//this->Approachlist.insert(0);

	}
	Vertex(double x, double y, double z)
	{
		this->id = -1;
		this->pos = S_Vec3d(x, y, z);
		//this->Approachlist.clear();
		//this->Approachlist.insert(0);

	}



};



//四阶矩阵 用于网格简化
typedef struct Matrix
{
	Matrix()
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				this->M[i][j] = 0;
			}
		}
	}
	double M[4][4];

}Matrix;


typedef struct Vec4
{
	Vec4()
	{
		for (int i = 0; i < 4; i++)
		{
			M[i]=0;
		}
	}
	Vec4(double _a, double _b, double _c, double _d)
	{
		M[0] = _a;
		M[1] = _b;
		M[2] = _c;
		M[3] = _d;
	}
	double M[4];

}Vec4;

//点组
class V_Group
{
public:
	V_Group(void);
	~V_Group(void);
	Vertex Vgroup[MAX_SIZE]; //千万注意这里点是[1,cur_max];
	int curnum;  //当前最大点id;
	bool isDeleted[MAX_SIZE];//判断某个顶点是否已经被删除
	int add_V(Vertex);//向group中增加一个顶点，返回顶点编号
	void del_V(int);//在group中删除id顶点
	int getCommonVertexNum(int, int);

};

//边堆
class E_Heap
{
public:
	E_Heap(void);
	~E_Heap(void);
	struct cmp {
		bool operator() (Edge X, Edge Y) {
			return X.deItav > Y.deItav;
		}
	};
	priority_queue<Edge, std::vector<Edge>, cmp> EdgeQueue;//边队列，按照误差从小到大排序
	vector<Edge> Edgelist;
	map<pair<int, int>, int> mapEdgeToID;//建立顶点到边的映射关系
	bool isDeleted[MAX_SIZE];//标记哪些边被删除

	int curnum;//边数量
	void addEdge(Edge&);//加边
	void delEdge(Edge);//删边
	Edge getMinD();//删除误差最小的边
};

//三维向量计算
S_Vec3d operator - (const S_Vec3d& A);
S_Vec3d operator + (const S_Vec3d& A, const S_Vec3d& B);
S_Vec3d operator - (const S_Vec3d& A, const S_Vec3d& B);
S_Vec3d operator * (const S_Vec3d& A, const S_Vec3d& B);
S_Vec3d operator / (const S_Vec3d& A, const S_Vec3d& B);
S_Vec3d operator + (const S_Vec3d& A, const double& B);
S_Vec3d operator - (const S_Vec3d& A, const double& B);
S_Vec3d operator * (const S_Vec3d& A, const double& B);
S_Vec3d operator / (const S_Vec3d& A, const double& B);

//四阶矩阵计算

Matrix operator + (const Matrix& A, const Matrix& B);

//点操作
//点的邻接表插入
void Ver_insert(Vertex* v1, int id);
//点的邻接表删除
void Ver_del(Vertex* v1, int id);
//是否邻接
bool Ver_isApproach(Vertex* v1, int id);

//点操作，判断点是否相邻
bool NearTo(Vertex* v1, Vertex* v2);
//三点计算面系数
Vec4 ComputeFace(Vertex* v1, Vertex* v2, Vertex* v3);


/*求边坍塌后点位置
*
*     q11 q12 q13 q14              0
*     q12 q22 q23 q24    *   v   = 0
*     q13 q23 q33 q34              0
*      0    0   0   1              1
*/

Vec4 gass(Matrix m, Vec4 Y);

struct E_cmp {
	bool operator() (Edge X, Edge Y) {
		return X.deItav > Y.deItav;
	}
};

typedef struct init_message
{
	int start;
	int end;
	int av_Del;
	int split;
	E_Heap *eheap;
	V_Group *vGroup;
	vector<Edge> idlist;
	priority_queue<Edge, std::vector<Edge>, E_cmp> EdgeQueue;
}init_message;

#endif // COMMON_H_INCLUDED
