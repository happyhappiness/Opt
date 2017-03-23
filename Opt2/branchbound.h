
#ifndef BRANCH_BOUND
#define BRANCH_BOUND

#include <vector>
#include <stack>
#include "simplex.h"

class BranchBound
{
public:
	BranchBound(int row, int col, std::ifstream &fin);

	bool solve(double &opt, std::vector<double> &vars);

private:
	int row;
	int col;

	std::stack<DualSimplex*> simplexs;
	double **matrix;

	double targetVal;
	std::vector<double> varval;
	std::vector<double> candicateVar;
	double lowBound;
	int noIntIndex;


	//���ص�һ���������±�  û���򷵻�-1
	int findNoIntIndex(std::vector<double> &vec);
};

#endif