
#ifndef BRANCH_BOUND
#define BRANCH_BOUND

#include <vector>
#include <stack>
#include "simplex.h"

class BranchBound
{
public:
	BranchBound(int row, int col, int num_edge, ElementType **matrix);

	bool solve(double low, double &opt, std::vector<ElementType> &vars);

private:
	int num_edge;
	int row;
	int col;

	std::stack<DualSimplex*> simplexs;
	ElementType **matrix;

	double targetVal;
	std::vector<ElementType> varval;
	std::vector<ElementType> candicateVar;
	double lowBound;
	int noIntIndex;


	//返回第一个非整数下标  没有则返回-1
	int findNoIntIndexFromBegin(std::vector<ElementType> &vec);
	int findNoIntIndexFromEnd(std::vector<ElementType> &vec);
	int findNoIntIndex(std::vector<ElementType> &vec);

};

#endif