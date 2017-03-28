
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
	DualSimplex(const DualSimplex &oldSimplex, bool& hasFound);

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

	// clear old matrix for update simplex
	void clearOldMatrix()
	{
		for (int i = 0; i < num_row; i++)
			delete[]matrix[i];
		delete[]matrix; 
	}

	void setMatrix(int row, int col, ElementType **data);
	void readMatrix(int row, int col, std::ifstream & file);
	void outputMatrix(std::ofstream & file);
	void outputMatrix();

	//  execute solve
	bool solveMinProblemWithDual(double &best, std::vector<ElementType>& variableValue);


	// add new constraint to the old simplex
	bool updateSimplex();

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


	// get float part of the non-integeral element
	ElementType getFloatPart(ElementType element);
	// find non-int element in b column(the first column)
	// prefer the one near to 0.5
	int getNonIntegralB();
	// calculate the new constraints of cutting plane
	void getConstraint(int rowIndex, std::vector<ElementType>& constraint);


	int num_variable;// variables
	int num_row;// row = constraints + objectives
	int num_col;// col = variables + relaxation + b
	ElementType **matrix; // matrix for data
	std::vector<int> bv;// stroe basic value

};

#endif