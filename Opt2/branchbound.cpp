
#include "branchbound.h"
#include <iostream>
using namespace std;

BranchBound::BranchBound(int row, int col, int num_edge, ElementType **matrix)
{
	this->row = row;
	this->col = col;
	this->num_edge = num_edge;
	this->matrix = matrix;
}

bool BranchBound::solve(double low, double &opt, std::vector<ElementType> &vars)
{
	lowBound = low;

	DualSimplex *simplex0 = new DualSimplex();
	simplex0->setMatrix(row, col, matrix);
	simplexs.push(simplex0);

	int problemCnt = 0;

	while (!simplexs.empty())
	{
		//cout << simplexs.size() << endl;
		DualSimplex *curSimplex = simplexs.top();
		simplexs.pop();

		problemCnt++;
		if (!curSimplex->solveMinProblemWithDual(targetVal, varval))
		{
			delete curSimplex;
			continue;
		}
		
		if (targetVal <= lowBound)
		{
			delete curSimplex;
			continue;
		}
		noIntIndex = findNoIntIndex(varval);
		if (noIntIndex == -1)
		{
			delete curSimplex;
			lowBound = targetVal;
			candicateVar = varval;
			cout << lowBound << endl;}
		
		else
		{
			double b = varval[noIntIndex];
			//左分支
			DualSimplex *leftSimplex = new DualSimplex(*curSimplex, noIntIndex, (int)b, true);
			//右分支
			DualSimplex *rightSimplex = new DualSimplex(*curSimplex, noIntIndex, (int)b + 1, false);
			
			simplexs.push(leftSimplex);
			simplexs.push(rightSimplex);
			
			delete curSimplex;
		}
	}
	cout <<"problemCnt"<< problemCnt << endl;
	if (candicateVar.size() != 0)
	{
		opt = lowBound;
		vars = candicateVar;
		return true;
	}
	else
	{
		return false;
	}
}


//从前往后找
int BranchBound::findNoIntIndexFromBegin(std::vector<ElementType> &vec)
{
	int total = vec.size();
	for (int i = num_edge; i < total; i++)
	{
		double val = vec[i];
		if (abs(round(val) - val) > 0.00001)
		{
			return i;
		}
	}
	return -1;
}

//从后往前找
int BranchBound::findNoIntIndexFromEnd(std::vector<ElementType> &vec)
{
	int total = vec.size();
	for (int i = total - 1; i >= num_edge; i--)
	{
		double val = vec[i];
		if (abs(round(val) - val) > 0.00001)
		{
			return i;
		}
	}
	return -1;
}

//寻找偏离整数最大的
int BranchBound::findNoIntIndex(std::vector<ElementType> &vec)
{
	int total = vec.size();
	int finded = -1;
	double away = 0.00001;
	for (int i = num_edge; i < total; i++)
	{
		double val = vec[i];
		double thisAway = abs(round(val) - val);
		if (thisAway > away)
		{
			away = thisAway;
			finded = i;
		}
	}
	return finded;
}