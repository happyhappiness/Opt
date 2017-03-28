#include "simplex.h"
#include <Windows.h>
using namespace std;


DualSimplex::DualSimplex(const DualSimplex &oldSimplex, int varIndex, int b, bool less)
{
	num_row = oldSimplex.num_row + 1;
	num_variable = oldSimplex.num_variable;
	num_col = oldSimplex.num_col + 1; // b + variable + relaxation

	//申请空间
	matrix = new ElementType*[num_row];
	for (int i = 0; i < num_row; i++)
	{
		matrix[i] = new ElementType[num_col];
	}

	//把父问题的矩阵复制过来
	for (int i = 0; i < oldSimplex.num_row; i++)
	{
		for (int j = 0; j < oldSimplex.num_col; j++)
		{
			matrix[i][j] = oldSimplex.matrix[i][j];
		}
	}

	//添加约束和松弛变量
	int lastRow = num_row - 1;
	int lastCol = num_col - 1;
	// last col is 0
	for (int i = 0; i < lastRow; i++)
	{
		matrix[i][lastCol] = 0;//relax
	}
	// last row is 0
	for (int i = 0; i < lastCol; i++)
	{
		matrix[lastRow][i] = 0;
	}
	matrix[lastRow][0] = less ? b : -b;
	matrix[lastRow][varIndex + 1] = less ? 1 : -1;
	matrix[lastRow][lastCol] = 1; //relax

	//把基变量编号复制过来
	bv = oldSimplex.bv;
	bv.push_back(num_col - 2);

	//基变量冲突吗？
	int totalBv = bv.size();
	for (int i = 0; i < totalBv; i++)
	{
		if (bv[i] == varIndex)
		{
			//有冲突  把约束的相应变量消除
			if (less)
			for (int j = 0; j < num_col; j++)
				matrix[lastRow][j] -= matrix[i + 1][j];
			else
			for (int j = 0; j < num_col; j++)
				matrix[lastRow][j] += matrix[i + 1][j];
			break;
		}
	}
}


void DualSimplex::setMatrix(int row, int col, ElementType **data)
{
	clear();
	num_row = row;
	num_variable = col - 1;

	num_col = 1 + num_variable + num_row - 1; //b + variable + relaxation
	//cout << "reading from data matrix: " << num_row << ", col:" << col << endl;

	matrix = new ElementType*[num_row];
	for (int i = 0; i < num_row; i++)
	{
		// read b and variables
		matrix[i] = new ElementType[num_col];
		for (int j = 0; j < col; j++)
		{
			matrix[i][j] = data[i][j];
		}
	}
	//设置基变量
	for (int i = num_variable; i < num_col - 1; i++)
	{
		bv.push_back(i);
	}
	// add relexation
	addRelaxVars();
}

void DualSimplex::readMatrix(int row, int col, std::ifstream & file)
{
	num_row = row;
	num_variable = col - 1;

	num_col = 1 + num_variable + num_row - 1; // b + variable + relaxation
	cout << "reading from file row: " << row << ", col:" << col << endl;

	matrix = new ElementType*[num_row];
	//matrix = (double**)malloc(num_col*num_row*sizeof(double));
	for (int i = 0; i < num_row; i++)
	{
		// read b and variables
		matrix[i] = new ElementType[num_col];
		for (int j = 0; j < col; j++)
		{
			file >> matrix[i][j];
		}
	}

	// add relexation
	addRelaxVars();
}

void DualSimplex::outputMatrix(std::ofstream & file)
{
	for (int i = 0; i < num_row; i++)
	{
		// read b and variables
		for (int j = 0; j < num_col; j++)
		{
			file << matrix[i][j] << "\t";
		}
		file << endl;
	}

}

void DualSimplex::outputMatrix()
{
	for (int i = 0; i < num_row; i++)
	{
		// read b and variables
		for (int j = 0; j < num_col; j++)
		{
			cout << matrix[i][j] << "\t";
		}
		cout << endl;
	}
}

// add relax variables for less constraint
void DualSimplex::addRelaxVars()
{
	// the first row for the relaxation is 0
	for (int j = 1; j < num_row; j++)
	{
		matrix[0][num_variable + j] = 0;
	}

	// num_relaxation + num_variable + 1 = num_col (num_row = num_relaxation + 1)
	// from the second row to last row: add relaxation (the first row is objective, all constraints is less)
	for (int i = 1; i < num_row; i++)
	{
		// from the column with index num_variable + 1, to num_col: add relaxation (the first colunm is b)
		for (int j = 1; j < num_row; j++)
		{
			matrix[i][num_variable + j] = (i == j);
		}
	}
}

// find the row to pivot (most negative row)
int DualSimplex::findPivotRow()
{
	// find the most negative one in b
	double lowest = 0;
	int pivot_row = -1;
	for (int i = 1; i < num_row; i++)
	{
		double val = matrix[i][0];
		//// accuracy
		//if (abs(round(val) - val) < 1e-8)
		//	matrix[i][0] = round(val);
		if (matrix[i][0] < lowest)
		{
			lowest = matrix[i][0];
			pivot_row = i;
		}
	}

	return pivot_row;
}

// find the column to pivot (ratio with least abs value)
int DualSimplex::findPivotCol(int row)
{
	double largest = -INFINITY;
	int pivot_col = -1;
	double now_ratio = 0;
	// the first col is b
	for (int j = 1; j < num_col; j++)
	{
		// try to make b to be positive, so choose the negative col && the first row must be non-negative
		if (matrix[row][j] >= 0)
			continue;

		now_ratio = matrix[0][j] / matrix[row][j];
		// choose the negative col with min abs ratio, max negative
		if (now_ratio > largest)
		{
			largest = now_ratio;
			pivot_col = j;
		}
	}
	return pivot_col;
}

// execute pivotion
void DualSimplex::pivot(int row, int col)
{
	bv[row - 1] = col - 1;
	double cofficiency = matrix[row][col];

	// cast the cofficiency of pivot into 1
	for (int j = 0; j < num_col; j++)
		matrix[row][j] /= cofficiency;

	// cast cofficiency of pivot of the remained row into 0
	double multiplier;
	for (int i = 0; i < num_row; i++) {
		if (i == row)
			continue;
		// update each row 
		multiplier = matrix[i][col];
		if (multiplier == 0)
			continue;

		for (int j = 0; j < num_col; j++)
		{
			matrix[i][j] -= multiplier * matrix[row][j];
		}
	}
}

// get optimization variables and value
bool DualSimplex::getOptimization(double &best, vector<ElementType>& variableValue)
{
	variableValue.clear();
	for (int i = 0; i < num_variable; i++)
	{
		variableValue.push_back(0);
	}
	int totalBv = bv.size();
	for (int i = 0; i < totalBv; i++)
	{
		int bvNum = bv[i];
		if (bvNum < num_variable)
			variableValue[bvNum] = matrix[i + 1][0];
	}
	best = matrix[0][0];

	return true;
}


//  execute solve (after readMatrix or setMatrix)
bool DualSimplex::solveMinProblemWithDual(double &best, vector<ElementType>& variableValue)
{
	/*LARGE_INTEGER freq;
	LARGE_INTEGER start_t, stop_t;
	double exe_time;
	QueryPerformanceFrequency(&freq);

	QueryPerformanceCounter(&start_t);*/

	// iteration at most loop times
	int pivotRow, pivotCol;
	bool flag = false;
	while (1)
	{
		//step1: find pivot row
		pivotRow = findPivotRow();
		// no negative b value, already the best solution
		if (pivotRow == -1)
		{
			// get optimized solution
			getOptimization(best, variableValue);
			flag = true;
			break;
		}

		pivotCol = findPivotCol(pivotRow);
		// exist some row that can never have positive b value
		if (pivotCol == -1)
		{
			//ofstream fout("test.txt");
			//outputMatrix(fout);
			//fout.close();

			flag = false;
			break;
		}

		//step2 do pivot
		pivot(pivotRow, pivotCol);
		
	}
	return flag;
}

// get float part of the non-integeral element
ElementType DualSimplex::getFloatPart(ElementType element)
{
	// check accuracy firstly
	return abs(element - round(element)) > accuracy ? (element - floor(element)) : 0;
}

// find non-int element in b column(the first column)
// prefer the one near to 0.5
int DualSimplex::getNonIntegralB()
{
	// find nearest float part
	ElementType nearest = 0.5;
	ElementType nowB, delta;
	int rowIndex = -1;
	// satrt from row with b value
	for (int i = 1; i < this->num_row; i++)
	{
		nowB = getFloatPart(this->matrix[i][0]);
		delta = abs(0.5 - nowB);
		if (delta < nearest)
		{
			nearest = delta;
			rowIndex = i;
		}
	}

	return rowIndex;
}
// calculate the new constraints of cutting plane
void DualSimplex::getConstraint(int rowIndex, vector<ElementType>& constraint)
{
	// b value(negative)
	ElementType floatB = getFloatPart(this->matrix[rowIndex][0]);
	constraint.push_back(-1 * floatB);

	ElementType floatPart;
	int sign;
	// calculate opssite sign * float part for each variable
	for (int j = 1; j < this->num_col; j++)
	{
		floatPart = getFloatPart(this->matrix[rowIndex][j]);
		if (floatPart == 0)
		{
			constraint.push_back(0);
			continue;
		}
		// calculate the opposite sign
		sign = this->matrix[rowIndex][j] > 0 ? -1 : 1;
		constraint.push_back(sign * floatPart);
	}

}

// add new constraint to dualSimplex by finding the non-integeral b value
DualSimplex::DualSimplex(const DualSimplex &oldSimplex, bool& hasFound)
{
	// get non integeral b value
	int rowIndex = this->getNonIntegralB();
	if (rowIndex == -1)
	{
		hasFound = false;
		return;
	}

	hasFound = true;
	// get the new constraint
	vector<ElementType> constraint;
	this->getConstraint(rowIndex, constraint);

	num_row = oldSimplex.num_row + 1;
	num_variable = oldSimplex.num_variable;
	num_col = oldSimplex.num_col + 1;
	int lastRow = num_row - 1;
	int lastCol = num_col - 1;

	// ask for space
	matrix = new ElementType*[num_row];
	for (int i = 0; i < num_row; i++)
	{
		matrix[i] = new ElementType[num_col];
	}
	// copy the old matrix
	for (int i = 0; i < lastRow; i++)
	{
		for (int j = 0; j < lastCol; j++)
		{
			matrix[i][j] = oldSimplex.matrix[i][j];
		}
	}

	// add the last row and corresponding relaxation
	// last col is 0
	for (int i = 0; i < lastRow; i++)
	{
		matrix[i][lastCol] = 0;//relax
	}
	// last row is decided by constraint
	for (int i = 0; i < lastCol; i++)
	{
		matrix[lastRow][i] = constraint[i];
	}
	matrix[lastRow][lastCol] = 1; //relax

	// copy basic variable
	bv = oldSimplex.bv;
	bv.push_back(lastCol - 1);

	////基变量冲突吗？impossible
	//int totalBv = bv.size();
	//for (int i = 0; i < totalBv; i++)
	//{
	//	if (bv[i] == varIndex)
	//	{
	//		//有冲突  把约束的相应变量消除
	//		if (less)
	//			for (int j = 0; j < num_col; j++)
	//				matrix[lastRow][j] -= matrix[i + 1][j];
	//		else
	//			for (int j = 0; j < num_col; j++)
	//				matrix[lastRow][j] += matrix[i + 1][j];
	//		break;
	//	}
	//}

}


bool DualSimplex::updateSimplex()
{
	// get non integeral b value
	int rowIndex = this->getNonIntegralB();
	if (rowIndex == -1)
	{
		return false;
	}

	// get the new constraint
	vector<ElementType> constraint;
	this->getConstraint(rowIndex, constraint);

	int newRow = num_row + 1;
	int newCol = num_col + 1;
	// get new matrix
	// ask for space
	ElementType ** tempMatrix = new ElementType*[newRow];
	for (int i = 0; i < newRow; i++)
	{
		tempMatrix[i] = new ElementType[newCol + 1];
	}
	// copy the old matrix
	for (int i = 0; i < num_row; i++)
	{
		for (int j = 0; j < num_col; j++)
		{
			tempMatrix[i][j] = matrix[i][j];
		}
	}
	// add the last row and corresponding relaxation
	// last col is 0
	for (int i = 0; i < num_row; i++)
	{
		tempMatrix[i][num_col] = 0;//relax
	}
	// last row is decided by constraint
	for (int i = 0; i < num_col; i++)
	{
		tempMatrix[num_row][i] = constraint[i];
	}
	tempMatrix[num_row][num_col] = 1; //relaxation

	// update basic variable
	bv.push_back(num_col - 1);
	// update the matrix
	clearOldMatrix();
	matrix = tempMatrix;
	// update row and col
	num_row = newRow;
	num_col = newCol;

	return true;
}

//// add new constraint to the old simplex
//bool DualSimplex::updateSimplex(DualSimplex &newSimplex, vector<ElementType>& constraint)
//{
//
//}