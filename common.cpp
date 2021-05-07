#include "common.h"

bool NearTo(Vertex* v1, Vertex* v2)
{
	int id2;
	id2 = v2->id;
	for (set<int>::iterator it1 = v1->Approachlist.begin(); it1 != v1->Approachlist.end(); it1++)
	{
		if (*it1 == id2)
		{
			return true;
		}
	}
	return false;
	//test
	
}

Vec4 ComputeFace(Vertex* v1, Vertex* v2, Vertex* v3)
{
	double x1, x2, x3;
	double y1, y2, y3;
	double z1, z2, z3;
	x1 = v1->pos.x;
	x2 = v2->pos.x;
	x3 = v3->pos.x;
	y1 = v1->pos.y;
	y2 = v2->pos.y;
	y3 = v3->pos.y;
	z1 = v1->pos.z;
	z2 = v2->pos.z;
	z3 = v3->pos.z;
	float A;
	float B;
	float C;
	float D;
	A = (y3 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
	B = (x3 - x1) * (z2 - z1) - (x2 - x1) * (z3 - z1);
	C = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
	D = -(A * x1 + B * y1 + C * z1);
	return Vec4(A,B,C,D);
}

Vec4 gass(Matrix m, Vec4 Y)
{
	Vec4 v;
	for (int i = 0; i < 4; i++) {
		int j = 0;
		while (j < 4 && fabs(m.M[i][j]) < EPS)
			j++;
		if (j == 4)
			continue;
		for (int k = 0; k < 4; k++) {
			if (k != i) {
				double rate = m.M[k][j] / m.M[i][j];
				for (int l = 0; l < 4; l++)
					m.M[k][l] -= m.M[i][l] * rate;
				v.M[k] -= v.M[i] * rate;
			}
		}
	}
	Vec4 X;
	for (int i = 0; i < 4; i++) {
		int j = 0;
		while (j < 4 && fabs(m.M[i][j]) < EPS)
			j++;
		if (j == 4)
			return Vec4(0, 0, 0, -1);
		X.M[i] = v.M[i] / m.M[i][j];
	}
	return X;




}

//
S_Vec3d operator-(const S_Vec3d& A)
{
	return S_Vec3d(-A.x,-A.y,-A.z);
}

S_Vec3d operator+(const S_Vec3d& A, const S_Vec3d& B)
{
	return S_Vec3d(A.x + B.x, A.y + B.y, A.z + B.z);
}

S_Vec3d operator-(const S_Vec3d& A, const S_Vec3d& B)
{
	return S_Vec3d(A.x - B.x, A.y - B.y, A.z - B.z);
}

S_Vec3d operator*(const S_Vec3d& A, const S_Vec3d& B)
{
	return S_Vec3d(A.x * B.x , A.y * B.y , A.z * B.z );
}

S_Vec3d operator/(const S_Vec3d& A, const S_Vec3d& B)
{
	return S_Vec3d(A.x / B.x, A.y / B.y, A.z / B.z);
}

S_Vec3d operator+(const S_Vec3d& A, const double& B)
{
	return S_Vec3d(A.x+B,A.y+B,A.z+B);
}

S_Vec3d operator-(const S_Vec3d& A, const double& B)
{
	return S_Vec3d(A.x - B, A.y - B, A.z - B);
}

S_Vec3d operator*(const S_Vec3d& A, const double& B)
{
	return S_Vec3d(A.x * B, A.y * B, A.z * B);
}

S_Vec3d operator/(const S_Vec3d& A, const double& B)
{
	return S_Vec3d(A.x / B, A.y / B, A.z / B);
}

Matrix operator+(const Matrix& A, const Matrix& B)
{
	Matrix m;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			m.M[i][j] = A.M[i][j] + B.M[i][j];
	return m;
}

//
void Ver_insert(Vertex* v1, int id)
{
	v1->Approachlist.insert(id);
}

void Ver_del(Vertex* v1, int id)
{
	v1->Approachlist.erase(id);
}

bool Ver_isApproach(Vertex* v1, int id)
{
	return(v1->Approachlist.count(id) > 0);
}


V_Group::V_Group(void)
{
	curnum = 0;
	for (int i = 0; i < MAX_SIZE; i++)
		isDeleted[i] = false;

}

V_Group::~V_Group(void)
{
}

int V_Group::add_V(Vertex v)
{
	curnum++;
	v.id = curnum;
	v.isbound=0;
	Vgroup[curnum] = v;
	return curnum;
}

void V_Group::del_V(int id)
{
	if (id >= MAX_SIZE)
	{
		return;
	}
	isDeleted[id] = true;
	for (set<int>::iterator it = Vgroup[id].Approachlist.begin(); it != Vgroup[id].Approachlist.end(); it++) {
		Ver_del(&Vgroup[(*it)], id);
	}
}

int V_Group::getCommonVertexNum(int u, int v)
{
	int cnt = 0;
	for (set<int>::iterator it = Vgroup[u].Approachlist.begin();
		it != Vgroup[u].Approachlist.end(); it++) {//在u的邻接点中遍历
		if (NearTo( Vgroup+v, Vgroup + *it)) {
			cnt++;
		}

	}
	return cnt;
}

E_Heap::E_Heap(void)
{
	curnum = 0;
	for (int i = 0; i < MAX_SIZE; i++)
		isDeleted[i] = false;

}

E_Heap::~E_Heap(void)
{
}

void E_Heap::addEdge(Edge &e)
{
	curnum++;
	e.id = curnum;
	int u = e.v1;
	int v = e.v2;
	mapEdgeToID[make_pair(u, v)] = curnum;
	e.isbound=0;
	Edgelist.push_back(e);
}

void E_Heap::delEdge(Edge e)
{
	int u = e.v1;
	int v = e.v2;
	int ID= mapEdgeToID[make_pair(u, v)];
	isDeleted[ID] = true;
}


Edge E_Heap::getMinD()
{
	if (EdgeQueue.size() <= 0) {
		return Edge(0, 0);
	}
	while (isDeleted[EdgeQueue.top().id]) {
		EdgeQueue.pop();
	}
	Edge e = EdgeQueue.top();
	EdgeQueue.pop();
	return e;
	return Edge();
}
