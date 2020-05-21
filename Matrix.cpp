#include <stdexcept>
#include <iostream>
#include "Matrix.hpp"
#include <cmath>
using namespace std;

void Matrix::init()
{
	if ((r == 0) && (c == 0))
	{
		valid = false;
		return;
	}

	for(int i = 0; i < r; i++) {
        matrix.push_back({});
	    for(int j = 0; j < c; j++) {
            matrix[i].push_back(0);
	    }
	}
}

void Matrix::clean()
{
	c = r = 0;
	valid = false;
}

Matrix::Matrix()
{
	r = 0;
	c = 0;
	valid = false;
}

Matrix::Matrix(size_t cols)
{
	r = 1;
	c = cols;
	valid = c != 0;

	init();
}

Matrix::Matrix(size_t rows, size_t cols)
{
	r = rows;
	c = cols;
	
	valid = r * c != 0;

	init();
}

Matrix::~Matrix()
{
	clean();
}

Matrix::Matrix(const Matrix& mat)
{
	if (!mat.valid)
	{
		return;
	}
	r = mat.r;
	c = mat.c;
	valid = mat.valid;
	init();

    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            matrix[i][j] = mat.matrix[i][j];
}

bool Matrix::isValid()
{
	return valid;
}

Matrix Matrix::operator*(const Matrix& mat) const
{
	if (c != mat.r)
	{
		Matrix res;
		res.valid = false;
		return res;
	}
	if ((r == 0) || (c == 0))
		return Matrix();

	Matrix cps(r, mat.c);
	for (int i = 0; i < r; i++)
		for (int j = 0; j < mat.c; j++)
		{
			cps.matrix[i][j] = 0;
			for (int k = 0; k < c; k++)
				cps.matrix[i][j] += matrix[i][k] * mat.matrix[k][j];
		}
	return cps;
}

Matrix Matrix::operator-(const Matrix& mat) const
{
	if ((r != mat.r) || (c != mat.c) || (r == 0) || (c == 0))
	{
		Matrix res;
		res.valid = false;
		return res;
	}

	Matrix summ(r, c);

	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			summ.matrix[i][j] = matrix[i][j] - mat.matrix[i][j];
	return summ;
}

Matrix Matrix::operator+(const Matrix& mat) const
{
	if ((r != mat.r) || (c != mat.c) || (r == 0) || (c == 0))
	{
		Matrix res;
		res.valid = false;
		return res;
	}

	Matrix summ(r, c);

	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			summ.matrix[i][j] = matrix[i][j] + mat.matrix[i][j];
	return summ;
}

Matrix Matrix::operator*(double value) const
{
	if((r == 0) || (c == 0))
		return Matrix();

	Matrix res(r, c);
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			res.matrix[i][j] = matrix[i][j] * value;
	return res;
}

Matrix Matrix::operator/(double value) const
{
	if ((value == 0) || (r == 0) || (c == 0))
		return Matrix();

	Matrix res(r, c);
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			res.matrix[i][j] = matrix[i][j] / value;
	return res;
}

Matrix& Matrix::operator=(const Matrix& mat)
{
	clean();
	r = mat.r;
	c = mat.c;
	valid = mat.valid;
	init();

	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			matrix[i][j] = mat.matrix[i][j];
	return *this;
}

Matrix& Matrix::operator*=(const Matrix& mat)
{
	if ((c != mat.r) || (r == 0) || (c == 0))
	{
		valid = false;
		return *this;
	}

	c = mat.c;

	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
		{
			double item = 0;
			for (int k = 0; k < c; k++)
				item += matrix[i][k] * mat.matrix[k][j];
			matrix[i][j] = item;
		}
	return *this;
}

Matrix& Matrix::operator+=(const Matrix& mat)
{
	if ((r != mat.r) || (c != mat.c) || (r == 0) || (c == 0))
	{
		valid = false;
		return *this;
	}

	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			matrix[i][j] += mat.matrix[i][j];

	return *this;
}

Matrix& Matrix::operator-=(const Matrix& mat)
{
	if ((r != mat.r) || (c != mat.c) || (r == 0) || (c == 0))
	{
		valid = false;
		return *this;
	}

	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			matrix[i][j] -= mat.matrix[i][j];

	return *this;
}

Matrix& Matrix::operator*=(double value)
{
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			matrix[i][j] *= value;

	return *this;
}

Matrix& Matrix::operator/=(double value)
{
	if (value == 0)
	{
		valid = false;
		return *this;
	}

	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			matrix[i][j] /= value;

	return *this;
}

Matrix operator*(double value, const Matrix& mat)
{
	Matrix res(mat);
	return res*value;
}

void Matrix::resize(size_t rows, size_t cols)
{
	if(valid)
		clean();
	r = rows;
	c = cols;
	init();
}

double& Matrix::coeffRef(size_t rowIdx, size_t colIdx)
{
	if ((rowIdx > r) || (colIdx > c))
		throw std::out_of_range("There is no such a cell in your matrix.");

	return matrix[rowIdx][colIdx];
}

const double& Matrix::coeffRef(size_t rowIdx, size_t colIdx) const
{
	if ((rowIdx > r) || (colIdx > c))
		throw std::out_of_range("message");

	return matrix[rowIdx][colIdx];
}

vector<vector<double>>& Matrix::data() {
	return matrix;
}

const vector<vector<double>>& Matrix::data() const {
	return matrix;
}

size_t Matrix::rows() const
{
	return (size_t)r;
}

size_t Matrix::cols() const
{
	return (size_t)c;
}

Matrix& Matrix::setIdentity()
{
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			if (i == j)
				matrix[i][j] = 1;
			else
				matrix[i][j] = 0;

	return *this;
}

Matrix& Matrix::setZero()
{
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			matrix[i][j] = 0;

	return *this;
}

Matrix& Matrix::setConstants(double value)
{
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			matrix[i][j] = value;

	return *this;
}

Matrix& Matrix::setIdentity(size_t rows, size_t cols)
{
	resize(rows, cols);
	return setIdentity();
}

Matrix& Matrix::setZero(size_t rows, size_t cols)
{
	resize(rows, cols);
	return setZero();
}

Matrix& Matrix::setConstants(size_t rows, size_t cols, double value)
{
	resize(rows, cols);
	return setConstants(value);
}

Matrix find_minor(const Matrix& pmatrix, int i, int j)
{
	int rows = pmatrix.rows();
	int cols = pmatrix.cols();

	Matrix minor(rows-1, cols-1);

	int b = 0;
	int c = 0;
	for (int k = 0; k < rows; k++)
	{
		if (k == i)
			continue;
		c = 0;
		for (int s = 0; s < cols; s++)
		{
			if (s == j)
				continue;
			minor.coeffRef(b, c) = pmatrix.coeffRef(k, s);
			c++;
		}
		b++;
	}

	return(minor);
}

Matrix Matrix::minor(int i, int j)
{
	Matrix minor(r - 1, c - 1);

	int b = 0;
	int l = 0;
	for (int k = 0; k < r; k++)
	{
		if (k == i)
			continue;
		l = 0;
		for (int s = 0; s < c; s++)
		{
			if (s == j)
				continue;
			minor.coeffRef(b, l) = matrix[k][s];
			l++;
		}
		b++;
	}

	return(minor);
}

bool Matrix::linaryDep(double* r1, double* r2, size_t size_row)
{
	double k = 0;
	for (int i = 0; i < size_row; i++)
	{
		if ((r1[i] != 0) && (k == 0))
			k = r2[i] / r1[i];

		if (r1[i] * k != r2[i])
			return false;
	}

	return true;
}

void Matrix::print()
{
	if ((r == 0) || (c == 0))
	{
		cout << "||" << endl;
		return;
	}

	cout.precision(2);

	cout << '\n';
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			double item = matrix[i][j];
			if (j == 0)
				cout << "|	" << item << "	";
			else if (j == (c - 1))
				cout << item << "	|";
			else
				cout << item << "	";
		}
		cout << '\n';
	}
	cout << '\n';
}

double* Matrix::solve()
{
	double* dets = nullptr;
	double* x = nullptr;

	Matrix system(r, c-1);
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c - 1; j++)
			system.matrix[i][j] = matrix[i][j];

	double det = system.det();

	if ((det == 0) || (r != (c - 1)))
		return x;

	dets	= new double[r];
	x		= new double[r];

	for (int i = 0; i < r; i++)
	{
		Matrix M_i(system);
		for (int j = 0; j < r; j++)
			M_i.matrix[j][i] = matrix[j][c - 1];
		double det_i = M_i.det();
		x[i] = det_i / det;
	}
	return x;
}

void Matrix::printSolution()
{
	double* s = solve();

	cout << endl;
	if (s == nullptr)
	{
		cout << "There is no one solution." << endl;
		return;
	}

	for (int i = 0; i < r; i++)
		cout << "X" << i + 1 << " = " << s[i] << endl;
	cout << endl;
}

Matrix Matrix::diag()
{
	Matrix diag_matrix(*this);	// Copy of our matrix (We must not crash our old matrix)

	int used_rows = 0;	// Iteratior for diag_matrix (the number of independent row)
	
	for (int j = 0; j < c; j++)
	{
		bool row = false;
		int i;
		for (i = used_rows; i < r; i++) // (***)
		{
			if (diag_matrix.matrix[i][j] != 0)
			{
				row = true;
				break;
			}
		}
		
		if (!row)
			continue;

		double first_item = diag_matrix.matrix[i][j];

		/*	We should select used rows. 
		 Therefore we replace [i] and [used_rows] rows 
		 and when we try to find new linery
		 independent row in (***) we do not consider
		 old rows by this technique.
		*/
		for (int s = 0; s < r; s++)
		{
			double item_i = diag_matrix.matrix[i][s];
			diag_matrix.matrix[i][s] = diag_matrix.matrix[used_rows][s];
			diag_matrix.matrix[used_rows][s] = item_i;
		}


		used_rows++;
		
		for (int s = used_rows; s < r; s++)
		{
			double koef = diag_matrix.matrix[s][j] / first_item;

			for (int k = 0; k < c; k++) 
				diag_matrix.matrix[s][k] -= diag_matrix.matrix[(used_rows - 1)][k] * koef;
		}
	}

	for (int i = used_rows+1; i < r; i++) // All other rows is linearly dependent.
		for (int j = 0; j < c; j++)
			diag_matrix.matrix[i][j] = 0;
	
	return diag_matrix;
}

Matrix Matrix::transpose()
{
	Matrix newm(c, r);
	
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++) 
			newm.coeffRef(j, i) = matrix[i][j];
	
	return newm;
}

double Matrix::det()
{
	double determinant = 0;

	if ((c != r) || (r == 0) || (c == 0)) {
        valid = false;
	    return (0);
    }

	if (r == 1)
		return matrix[0][0];
	else
	{
		for (int j = 0; j < c; j++) {
			Matrix minor(find_minor(*this, 0, j));
			determinant += pow((-1), j) * matrix[0][j] * minor.det();
		}
	}

	return(determinant);
}

Matrix Matrix::inverse()
{
	double k = det();

	if ((r != c) || (k == 0))
	{
		Matrix res;
		res.valid = false;
		return res;
	}
	if((r == 0) || (c == 0))
		return(Matrix());

	Matrix result = Matrix(*this);

	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
		{
			Matrix minor(find_minor(*this, j, i));
			result.matrix[i][j] = pow((-1), i + j) * minor.det() / k;
		}

	return(result);
}

Matrix Matrix::identity(size_t rows, size_t cols)
{
	return Matrix(rows, cols).setIdentity();
}

Matrix Matrix::zeros(size_t rows, size_t cols)
{
	return Matrix(rows, cols).setZero();
}

Matrix Matrix::constants(size_t rows, size_t cols, double value)
{
	return Matrix(rows, cols).setConstants(value);
}