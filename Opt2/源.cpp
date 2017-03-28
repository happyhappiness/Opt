
#include<fstream>
#include<iostream>
#include<vector>
#include<unordered_map>
#include<ctime>
#include<iomanip>

#include "simplex.h"
#include "branchbound.h"

using namespace std;

#define V 1500
#define E 11000
#define inf 99999999

struct Edge{
	int u, v, c, cost;
}edge[E];

vector<int> edgeIn[V];
vector<int> edgeOut[V];

int cnt = 0;
void addCenterEdge(int u, int v, int c, int cost)
{
	edge[cnt].u = u; edge[cnt].v = v; edge[cnt].cost = cost; edge[cnt].c = c;
	edgeIn[v].push_back(cnt);
	edgeOut[u].push_back(cnt);
	cnt++;

	edge[cnt].u = v; edge[cnt].v = u; edge[cnt].cost = cost; edge[cnt].c = c;
	edgeIn[u].push_back(cnt);
	edgeOut[v].push_back(cnt);
	cnt++;
}

unordered_map<int, int> neighborToEdge;
void addConsumerEdge(int u, int demand)
{
	edge[cnt].c = demand;
	neighborToEdge.insert(pair<int, int>(u, cnt));
	cnt++;
}


void setArray(ElementType array[], vector<int> &vec, int var)
{
	for (vector<int>::iterator i = vec.begin(); i != vec.end(); i++)
	{
		array[*i + 1] = var;
	}
}


void makeCons(char *path, int &totalRow, int &totalCol, int &num_edge, int &lowbound, ElementType **&matrix)
{
	ifstream fin(path);

	int vCnt, eCnt, cCnt;
	fin >> vCnt >> eCnt >> cCnt;

	totalRow = 1 + 2 * eCnt + 4 * vCnt;
	totalCol = 1 + 2 * eCnt + vCnt;
	num_edge = 2 * eCnt;
	matrix = new ElementType*[totalRow];
	for (int i = 0; i < totalRow; i++)
	{
		matrix[i] = new ElementType[totalCol];
	}

	int sCost;
	fin >> sCost;

	lowbound = sCost * cCnt;

	//内部网络
	int from, to, cap, cost;
	for (int i = 0; i < eCnt; i++)
	{
		fin >> from >> to >> cap >> cost;
		addCenterEdge(from, to, cap, cost);
	}

	//边缘网络
	int cId, vId, demand;
	int totalDemand = 0;
	for (int i = 0; i < cCnt; i++)
	{
		fin >> cId >> vId >> demand;
		totalDemand += demand;
		addConsumerEdge(vId, demand);
	}

	//*****************************构造约束********************************
	//首先是目标函数系数  前2*eCnt为边流量f的系数 后vCnt为节点是否为服务器s的系数
	int matrixRow = 0;
	matrix[matrixRow][0] = 0;
	for (int i = 0; i < 2 * eCnt; i++)
	{
		matrix[matrixRow][i+1] = edge[i].cost;
	}
	for (int i = 1; i <= vCnt; i++)
	{
		matrix[matrixRow][i + 2 * eCnt] = sCost;
	}

	int M = totalDemand;
	//容量约束
	matrixRow++;
	for (int i = 0; i < 2 * eCnt; i++)
	{
		memset(matrix[matrixRow], 0, sizeof(ElementType) * totalCol);
		matrix[matrixRow][0] = edge[i].c;
		matrix[matrixRow][i+1] = 1;
		matrixRow++;
	}


	//顶点0 1约束
	for (int i = 0; i < vCnt; i++)
	{
		memset(matrix[matrixRow], 0, sizeof(ElementType)* totalCol);
		matrix[matrixRow][0] = 1;
		matrix[matrixRow][i + 2 * eCnt + 1] = 1;
		matrixRow++;
	}


	//顶点流量平衡
	for (int i = 0; i < vCnt; i++)
	{
		memset(matrix[matrixRow], 0, sizeof(ElementType)* totalCol);
		matrix[matrixRow][0] = M;
		setArray(matrix[matrixRow], edgeIn[i], 1);
		matrix[matrixRow][2 * eCnt + i + 1] = M;
		matrixRow++;

		memset(matrix[matrixRow], 0, sizeof(ElementType)* totalCol);
		unordered_map<int, int>::iterator ite = neighborToEdge.find(i);
		int b = ite != neighborToEdge.end() ? edge[ite->second].c : 0;
		matrix[matrixRow][0] = b;
		setArray(matrix[matrixRow], edgeIn[i], 1);
		setArray(matrix[matrixRow], edgeOut[i], -1);
		matrix[matrixRow][2 * eCnt + i + 1] = -M;
		matrixRow++;

		memset(matrix[matrixRow], 0, sizeof(ElementType)* totalCol);
		matrix[matrixRow][0] = -b;
		setArray(matrix[matrixRow], edgeOut[i], 1);
		setArray(matrix[matrixRow], edgeIn[i], -1);
		matrix[matrixRow][2 * eCnt + i + 1] = -M;
		matrixRow++;
	}


	/*ofstream fout("opt.txt");
	for (int i = 0; i < totalRow; i++)
	{
		for (int j = 0; j < totalCol; j++)
		{
			fout << matrix[i][j] << " ";
		}
		fout << endl;
	}*/
	
}

int main()
{
	int row, col, num_edge, lowbound;
	ElementType **matrix = NULL;
	makeCons("case0.txt", row, col, num_edge, lowbound, matrix);

	//算法开始
	BranchBound branchBound(row, col, num_edge, matrix);
	double optimized;
	vector<ElementType> variables;
	if (!branchBound.solve(-lowbound-1, optimized, variables))
	{
		cout << "can't find best" << endl;
		return 0;
	}
	cout << "optimization value is " << optimized << endl;
	cout << "variables is " << endl;
	for (int i = 0; i < variables.size(); i++)
	{
		cout << setw(15) << variables[i];
	}
	cout << endl;

	return 0;
}

