#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

/*
*     �������ݽṹ������lod��mesh simplify
*     ��ά���� S_Vec3d /�� Edge /�� Vertex /��ά���� Matrix
*     ��group  V_Group  /   �߶�  E_Heap
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
//��ά����
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

//��
typedef struct Edge
{
	//ע�ⶥ����v1<v2
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
		this->deItav = 1000000; // ��һ�������ֵ
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
		this->deItav = 1000000; // ��һ�������ֵ
	}
	int id;//�߱��
	int v1, v2; //������
	S_Vec3d simp_v; //��̮�����
	double deItav; //��̮������
	int isbound;

}Edge;
//��
class Vertex
{
public:
	int id;
	S_Vec3d pos; //��ά����
	set<int> Approachlist; //�ڽӵ��
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



//�Ľ׾��� ���������
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

//����
class V_Group
{
public:
	V_Group(void);
	~V_Group(void);
	Vertex Vgroup[MAX_SIZE]; //ǧ��ע���������[1,cur_max];
	int curnum;  //��ǰ����id;
	bool isDeleted[MAX_SIZE];//�ж�ĳ�������Ƿ��Ѿ���ɾ��
	int add_V(Vertex);//��group������һ�����㣬���ض�����
	void del_V(int);//��group��ɾ��id����
	int getCommonVertexNum(int, int);

};

//�߶�
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
	priority_queue<Edge, std::vector<Edge>, cmp> EdgeQueue;//�߶��У���������С��������
	vector<Edge> Edgelist;
	map<pair<int, int>, int> mapEdgeToID;//�������㵽�ߵ�ӳ���ϵ
	bool isDeleted[MAX_SIZE];//�����Щ�߱�ɾ��

	int curnum;//������
	void addEdge(Edge&);//�ӱ�
	void delEdge(Edge);//ɾ��
	Edge getMinD();//ɾ�������С�ı�
};

//��ά��������
S_Vec3d operator - (const S_Vec3d& A);
S_Vec3d operator + (const S_Vec3d& A, const S_Vec3d& B);
S_Vec3d operator - (const S_Vec3d& A, const S_Vec3d& B);
S_Vec3d operator * (const S_Vec3d& A, const S_Vec3d& B);
S_Vec3d operator / (const S_Vec3d& A, const S_Vec3d& B);
S_Vec3d operator + (const S_Vec3d& A, const double& B);
S_Vec3d operator - (const S_Vec3d& A, const double& B);
S_Vec3d operator * (const S_Vec3d& A, const double& B);
S_Vec3d operator / (const S_Vec3d& A, const double& B);

//�Ľ׾������

Matrix operator + (const Matrix& A, const Matrix& B);

//�����
//����ڽӱ����
void Ver_insert(Vertex* v1, int id);
//����ڽӱ�ɾ��
void Ver_del(Vertex* v1, int id);
//�Ƿ��ڽ�
bool Ver_isApproach(Vertex* v1, int id);

//��������жϵ��Ƿ�����
bool NearTo(Vertex* v1, Vertex* v2);
//���������ϵ��
Vec4 ComputeFace(Vertex* v1, Vertex* v2, Vertex* v3);


/*���̮�����λ��
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
