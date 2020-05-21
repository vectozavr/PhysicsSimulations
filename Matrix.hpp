#pragma once

#include <vector>

class Matrix final {
private:
	std::vector<std::vector<double>> matrix;
	size_t r = 0;
	size_t c = 0;
	void init();
	void clean();
	bool valid = false;
public:
	Matrix();
	Matrix(size_t cols);
	Matrix(size_t rows, size_t cols);
	~Matrix();

	Matrix(const Matrix& mat);

	Matrix operator*(const Matrix& mat) const;
	Matrix operator-(const Matrix& mat) const;
	Matrix operator+(const Matrix& mat) const;

	Matrix operator*(double value) const;
	Matrix operator/(double value) const;

	Matrix& operator=(const Matrix& mat);
	Matrix& operator*=(const Matrix& mat);
	Matrix& operator+=(const Matrix& mat);
	Matrix& operator-=(const Matrix& mat);

	Matrix& operator*=(double value);
	Matrix& operator/=(double value);

	bool isValid();

	void resize(size_t rows, size_t cols);

	[[nodiscard]] const double& coeffRef(size_t rowIdx, size_t colIdx) const;
	double& coeffRef(size_t rowIdx, size_t colIdx);

    [[nodiscard]] const std::vector<std::vector<double>> & data() const;
    [[nodiscard]] std::vector<std::vector<double>> & data();

	size_t rows() const;
	size_t cols() const;

	Matrix& setIdentity();
	Matrix& setZero();
	Matrix& setConstants(double value);

	Matrix& setIdentity(size_t rows, size_t cols);
	Matrix& setZero(size_t rows, size_t cols);
	Matrix& setConstants(size_t rows, size_t cols, double value);

	Matrix diag();
	Matrix transpose();
	Matrix inverse();
	double det();

	//This is new useful functions:
	void print();			// Print matrix.
	void printSolution();	// Print solution for matrix (only for N*(N+1) sizes).
	double* solve();		// Solve matrix (only for N*(N+1) sizes).
	bool static linaryDep(double* r1, double* r2, size_t size_row); // Are r1 and r2 linery dependent.
	Matrix minor(int i, int j);
	//

	Matrix static identity(size_t rows, size_t cols);
	Matrix static zeros(size_t rows, size_t cols);
	Matrix static constants(size_t rows, size_t cols, double value);

	friend Matrix operator*(double value, const Matrix& mat);
};