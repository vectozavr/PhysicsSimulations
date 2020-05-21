#include "../googletest/googletest/include/gtest/gtest.h"
#include "../Matrix.hpp"

TEST(create, zero) {
	Matrix n(0, 0);
	EXPECT_EQ(n.cols(), 0);
	EXPECT_EQ(n.rows(), 0);
	EXPECT_ANY_THROW(n.coeffRef(0, 1));
	EXPECT_ANY_THROW(n.coeffRef(74, 85));
}

TEST(operators, all_operations_zero) {
	Matrix m1(0, 0);
	Matrix m2(0, 0);

	Matrix copy = m1;
	Matrix summ = m1 + m2;
	Matrix mult = m1 * m2;
	Matrix diff = m1 - m2;
	Matrix times = m1 * 10;
	Matrix div = m1 / 10;

	EXPECT_EQ(copy.cols(), 0);
	EXPECT_EQ(copy.rows(), 0);

	EXPECT_EQ(summ.cols(), 0);
	EXPECT_EQ(summ.rows(), 0);

	EXPECT_EQ(mult.cols(), 0);
	EXPECT_EQ(mult.rows(), 0);

	EXPECT_EQ(diff.cols(), 0);
	EXPECT_EQ(diff.rows(), 0);

	EXPECT_EQ(times.cols(), 0);
	EXPECT_EQ(times.rows(), 0);

	EXPECT_EQ(div.cols(), 0);
	EXPECT_EQ(div.rows(), 0);
}

TEST(operators, copy2x3) {
	Matrix m1(2, 3);

	m1.coeffRef(0, 0) = 1;
	m1.coeffRef(0, 1) = 0;
	m1.coeffRef(0, 2) = 9;
	m1.coeffRef(1, 0) = 2;
	m1.coeffRef(1, 1) = 3;
	m1.coeffRef(1, 2) = 7;

	Matrix copy = m1;

	EXPECT_EQ(copy.rows(), 2);
	EXPECT_EQ(copy.cols(), 3);

	EXPECT_EQ(copy.coeffRef(0, 0), 1);
	EXPECT_EQ(copy.coeffRef(0, 1), 0);
	EXPECT_EQ(copy.coeffRef(0, 2), 9);
	EXPECT_EQ(copy.coeffRef(1, 0), 2);
	EXPECT_EQ(copy.coeffRef(1, 1), 3);
	EXPECT_EQ(copy.coeffRef(1, 2), 7);
}

TEST(operators, summ2x3) {
	Matrix m1(2, 3);
	Matrix m2(2, 3);

	m1.coeffRef(0, 0) = 1;
	m1.coeffRef(0, 1) = 0;
	m1.coeffRef(0, 2) = 9;
	m1.coeffRef(1, 0) = 2;
	m1.coeffRef(1, 1) = 3;
	m1.coeffRef(1, 2) = 7;

	m2.coeffRef(0, 0) = 2;
	m2.coeffRef(0, 1) = 5;
	m2.coeffRef(0, 2) = 1;
	m2.coeffRef(1, 0) = 0;
	m2.coeffRef(1, 1) = 6;
	m2.coeffRef(1, 2) = 4;

	Matrix summ = m1 + m2;

	EXPECT_EQ(summ.rows(), 2);
	EXPECT_EQ(summ.cols(), 3);

	EXPECT_EQ(summ.coeffRef(0, 0), 3);
	EXPECT_EQ(summ.coeffRef(0, 1), 5);
	EXPECT_EQ(summ.coeffRef(0, 2), 10);
	EXPECT_EQ(summ.coeffRef(1, 0), 2);
	EXPECT_EQ(summ.coeffRef(1, 1), 9);
	EXPECT_EQ(summ.coeffRef(1, 2), 11);
}

TEST(operators, mult2x3) {
	Matrix m1(2, 3);
	Matrix m2(3, 2);
	Matrix m3(2, 3);

	m1.coeffRef(0, 0) = 1;
	m1.coeffRef(0, 1) = 0;
	m1.coeffRef(0, 2) = 9;
	m1.coeffRef(1, 0) = 2;
	m1.coeffRef(1, 1) = 3;
	m1.coeffRef(1, 2) = 7;

	m2.coeffRef(0, 0) = 3;
	m2.coeffRef(0, 1) = 6;
	m2.coeffRef(1, 0) = 8;
	m2.coeffRef(1, 1) = 1;
	m2.coeffRef(2, 0) = 2;
	m2.coeffRef(2, 1) = 4;

	m3.coeffRef(0, 0) = 1;
	m3.coeffRef(0, 1) = 0;
	m3.coeffRef(0, 2) = 9;
	m3.coeffRef(1, 0) = 2;
	m3.coeffRef(1, 1) = 3;
	m3.coeffRef(1, 2) = 7;

	Matrix mult = m1 * m2;
	Matrix invalid_mult = m1 * m3;
	
	EXPECT_EQ(mult.rows(), 2);
	EXPECT_EQ(mult.cols(), 2);
	EXPECT_FALSE(invalid_mult.isValid());

	EXPECT_EQ(mult.coeffRef(0, 0), 21);
	EXPECT_EQ(mult.coeffRef(0, 1), 42);
	EXPECT_EQ(mult.coeffRef(1, 0), 44);
	EXPECT_EQ(mult.coeffRef(1, 1), 43);
}

TEST(operators, diff2x3) {
	Matrix m1(2, 3);
	Matrix m2(2, 3);

	m1.coeffRef(0, 0) = 1;
	m1.coeffRef(0, 1) = 0;
	m1.coeffRef(0, 2) = 9;
	m1.coeffRef(1, 0) = 2;
	m1.coeffRef(1, 1) = 3;
	m1.coeffRef(1, 2) = 7;

	m2.coeffRef(0, 0) = 2;
	m2.coeffRef(0, 1) = 5;
	m2.coeffRef(0, 2) = 1;
	m2.coeffRef(1, 0) = 0;
	m2.coeffRef(1, 1) = 6;
	m2.coeffRef(1, 2) = 4;

	Matrix diff = m1 - m2;

	EXPECT_EQ(diff.rows(), 2);
	EXPECT_EQ(diff.cols(), 3);

	EXPECT_EQ(diff.coeffRef(0, 0), -1);
	EXPECT_EQ(diff.coeffRef(0, 1), -5);
	EXPECT_EQ(diff.coeffRef(0, 2), 8);
	EXPECT_EQ(diff.coeffRef(1, 0), 2);
	EXPECT_EQ(diff.coeffRef(1, 1), -3);
	EXPECT_EQ(diff.coeffRef(1, 2), 3);
}	

TEST(operators, times2x3) {
	Matrix m1(2, 3);

	m1.coeffRef(0, 0) = 1;
	m1.coeffRef(0, 1) = 0;
	m1.coeffRef(0, 2) = 9;
	m1.coeffRef(1, 0) = 2;
	m1.coeffRef(1, 1) = 3;
	m1.coeffRef(1, 2) = 7;

	Matrix times = m1 * 10;

	EXPECT_EQ(times.rows(), 2);
	EXPECT_EQ(times.cols(), 3);

	EXPECT_EQ(times.coeffRef(0, 0), 10);
	EXPECT_EQ(times.coeffRef(0, 1), 0);
	EXPECT_EQ(times.coeffRef(0, 2), 90);
	EXPECT_EQ(times.coeffRef(1, 0), 20);
	EXPECT_EQ(times.coeffRef(1, 1), 30);
	EXPECT_EQ(times.coeffRef(1, 2), 70);
}

TEST(operators, div2x3) {
	Matrix m1(2, 3);

	m1.coeffRef(0, 0) = 1;
	m1.coeffRef(0, 1) = 0;
	m1.coeffRef(0, 2) = 9;
	m1.coeffRef(1, 0) = 2;
	m1.coeffRef(1, 1) = 3;
	m1.coeffRef(1, 2) = 7;

	Matrix div = m1 / 10;

	EXPECT_EQ(div.rows(), 2);
	EXPECT_EQ(div.cols(), 3);

	EXPECT_EQ(div.coeffRef(0, 0), 0.1);
	EXPECT_EQ(div.coeffRef(0, 1), 0);
	EXPECT_EQ(div.coeffRef(0, 2), 0.9);
	EXPECT_EQ(div.coeffRef(1, 0), 0.2);
	EXPECT_EQ(div.coeffRef(1, 1), 0.3);
	EXPECT_EQ(div.coeffRef(1, 2), 0.7);
}

TEST(resize, zero) {
	Matrix n(0, 0);
	EXPECT_EQ(n.rows(), 0);
	EXPECT_EQ(n.cols(), 0);

	n.resize(10, 10);

	EXPECT_EQ(n.rows(), 10);
	EXPECT_EQ(n.cols(), 10);
}

TEST(resize, 2x5to5x5) {
	Matrix n(2, 5);
	EXPECT_EQ(n.rows(), 2);
	EXPECT_EQ(n.cols(), 5);

	n.resize(5, 5);

	EXPECT_EQ(n.rows(), 5);
	EXPECT_EQ(n.cols(), 5);
}

TEST(set, setIdentity) {
	Matrix n(7, 7);
	n.setIdentity();

	for (int i = 0; i < 7; i++)
		for (int j = 0; j < 7; j++)
			if (i == j)
				EXPECT_EQ(n.coeffRef(i, j), 1);
			else
				EXPECT_EQ(n.coeffRef(i, j), 0);
}

TEST(set, setZero) {
	Matrix n(7, 7);
	n.setZero();

	for (int i = 0; i < 7; i++)
		for (int j = 0; j < 7; j++)
			EXPECT_EQ(n.coeffRef(i, j), 0);
}

TEST(set, setConstant) {
	Matrix n(7, 7);
	n.setConstants(103);

	for (int i = 0; i < 7; i++)
		for (int j = 0; j < 7; j++)
			EXPECT_EQ(n.coeffRef(i, j), 103);
}

TEST(transpose, zero) {
	Matrix n(0, 0);
	Matrix res = n.transpose();

	EXPECT_EQ(res.rows(), 0);
	EXPECT_EQ(res.cols(), 0);
}

TEST(transpose, 7x7) {
	Matrix n(7, 7);
	Matrix res = n.setIdentity().transpose();
	
	for (int i = 0; i < 7; i++)
		for (int j = 0; j < 7; j++)
			EXPECT_EQ(res.coeffRef(i, j), n.coeffRef(i, j));
}

TEST(inverse, zero) {
	Matrix n(0, 0);
	Matrix res = n.inverse();

	EXPECT_EQ(res.rows(), 0);
	EXPECT_EQ(res.cols(), 0);
}

TEST(inverse, 5x5) {
	Matrix n(5, 5);
	n.setIdentity();
	n *= 2;
	Matrix inv(n.inverse());
	Matrix mult = inv * n;
	
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			if(i == j)
				EXPECT_EQ(mult.coeffRef(i, j), 1);
			else
				EXPECT_EQ(mult.coeffRef(i, j), 0);
}

TEST(det, 4x4) {
	Matrix n(4, 4);
	n.setIdentity();
	n *= 2;

	EXPECT_EQ(n.det(), 16);
}

TEST(det, 3x3) {
	Matrix n(3, 3);
	n.setIdentity();
	n *= 2;
	n.coeffRef(0, 1) = 1;
	n.coeffRef(0, 2) = 3;
	n.coeffRef(1, 2) = 3;

	EXPECT_EQ(n.det(), 8);
}

TEST(det, 1x1) {
	Matrix n(1, 1);
	n.coeffRef(0, 0) = 5;

	EXPECT_EQ(n.det(), 5);
}

TEST(diag, zero) {
	Matrix n(0, 0);
	Matrix d(n.diag());

	EXPECT_EQ(d.rows(), 0);
	EXPECT_EQ(d.cols(), 0);
}

TEST(diag, 3x3) {
	Matrix n(3, 3);
	n.coeffRef(0, 0) = 1;
	n.coeffRef(0, 1) = 2;
	n.coeffRef(0, 2) = 3;

	n.coeffRef(1, 0) = 5;
	n.coeffRef(1, 1) = 0;
	n.coeffRef(1, 2) = 1;

	n.coeffRef(2, 0) = 0;
	n.coeffRef(2, 1) = 0;
	n.coeffRef(2, 2) = 2;

	Matrix d = n.diag();

	EXPECT_EQ(d.coeffRef(0, 0), 1);
	EXPECT_EQ(d.coeffRef(0, 1), 2);
	EXPECT_EQ(d.coeffRef(0, 2), 3);
	EXPECT_EQ(d.coeffRef(1, 0), 0);
	EXPECT_EQ(d.coeffRef(1, 1), -10);
	EXPECT_EQ(d.coeffRef(1, 2), -14);
	EXPECT_EQ(d.coeffRef(2, 0), 0);
	EXPECT_EQ(d.coeffRef(2, 1), 0);
	EXPECT_EQ(d.coeffRef(2, 2), 2);
}

TEST(solve, 3x3) {
	Matrix n(3, 4);
	n.coeffRef(0, 0) = 1;
	n.coeffRef(0, 1) = 2;
	n.coeffRef(0, 2) = 3;

	n.coeffRef(1, 0) = 5;
	n.coeffRef(1, 1) = 0;
	n.coeffRef(1, 2) = 1;

	n.coeffRef(2, 0) = 0;
	n.coeffRef(2, 1) = 0;
	n.coeffRef(2, 2) = 2;

	n.coeffRef(0, 3) = 1;
	n.coeffRef(1, 3) = 1;
	n.coeffRef(2, 3) = 1;

	double* s = n.solve();
	EXPECT_EQ(s[0], 0.1);
	EXPECT_EQ(s[1], -0.3);
	EXPECT_EQ(s[2], 0.5);
}

TEST(mixComplexTest, 7x7)
{
	Matrix matr(7, 7);

	for (int i = 0; i < 7; i++)
		for (int j = 0; j < 7; j++)
			matr.coeffRef(i, j) = rand();

	Matrix inv = matr.inverse();
	Matrix comp = inv * matr;

	for (int i = 0; i < 7; i++)
		for (int j = 0; j < 7; j++)
		{
			if (i == j)
				EXPECT_TRUE(abs(comp.coeffRef(i, j) - 1) < 0.000001);
			else
				EXPECT_TRUE(abs(comp.coeffRef(i, j)	   ) < 0.000001);
		}
}

// Liza' tests

// Liza' tests

int main(int argc, char *argv[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}