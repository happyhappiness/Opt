
#ifndef DUAL_SIMPLEX
#define DUAL_SIMPLEX

#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include "util.h"


class DualSimplex
{
public:
	DualSimplex(){}
	DualSimplex(const DualSimplex &oldSimplex, int varIndex, int b, bool less);


	~DualSimplex()
	{
		clear();
	}

	void clear()
	{
		for (int i = 0; i < num_row; i++)
			delete[]matrix[i];
		delete[]matrix;
		bv.clear();
		num_variable = 0;
		num_row = 0;
		num_col = 0;
	}


	void setMatrix(int row, int col, ElementType **data);
	void readMatrix(int row, int col, std::ifstream & file);
	void outputMatrix(std::ofstream & file);
	void outputMatrix();

	//  execute solve
	bool solveMinProblemWithDual(double &best, std::vector<ElementType>& variableValue);
	

private:
	// add relax variables for less constraint
	void addRelaxVars();
	// find the row to pivot (most negative row)
	int findPivotRow();
	// find the column to pivot (ratio with least abs value)
	int findPivotCol(int row);
	// execute pivotion
	void pivot(int row, int col);
	// get optimization variables and value
	bool getOptimization(double &best, std::vector<ElementType>& variableValue);

	int num_variable;// variables
	int num_row;// row = constraints + objectives
	int num_col;// col = variables + relaxation + b
	ElementType **matrix; // matrix for data
	std::vector<int> bv;

	bool infeasible;
};

#endif