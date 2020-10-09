//
// Created by Иван Ильин on 08.10.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

int main() {

    // MOTION OF PLANETS
    ComplexPlot y_A;
    ComplexPlot y_SW;
    ComplexPlot y_SW_err;

    // range of integration
    double x_0 = 0;
    double x_N = PI;
    // initial conditions
    double y_0 = 1;
    double y_N = 0.5f;
    // integration step
    double h = 0.001;

    // the number of steps (calculated automatically)
    int N = (x_N - x_0) / h;

    vector<vector<double>> matrix;
    vector<double> d(N);
    vector<double> y(N);
    // Ay = d, here we fullfit the A matrix and d vector
    for(int i = 0; i < N; i++) {
        matrix.push_back({});
        for(int j = 0; j < N; j++) {
            if(j == i-1) {
                matrix[i].push_back(1.0f);
            } else if(j == i) {
                matrix[i].push_back(-2.0f);
            } else if(j == i + 1) {
                matrix[i].push_back(1.0f);
            } else {
                matrix[i].push_back(0.0f);
            }
        }
        d[i] = (h*h*sin(h*i));
    }

    // taking into account the initial conditions in A matrix
    matrix[0][0] = 1;
    matrix[0][1] = 0;

    matrix[N-1][N-1] = 1;
    matrix[N-1][N-2] = 0;
    // taking into account the initial conditions in d vector
    d[0] = y_0;
    d[N-1] = y_N;
    // back propagation
    for(int i = 1; i < N; i++) {
        double ksi = matrix[i][i-1] / matrix[i-1][i-1];
        matrix[i][i-1] = 0;
        matrix[i][i] -= ksi * matrix[i-1][i];
        d[i] -= ksi * d[i-1];
    }

    y[N-1] = d[N-1] / matrix[N-1][N-1];
    for(int i = N-2; i >= 0; i--)
        y[i] = (d[i] - matrix[i][i+1] * y[i+1]) / matrix[i][i];
    // save y in .dat file
    for(int i = 0; i < N; i++)
        y_SW.push(h*i, y[i]);

    saveVectorPoint2DToFile(y_SW.real(), "y_SW.dat");

    // analytical solution
    for(int i = 0; i < N; i++) {
        y_A.push(h*i, -sin(h*i) + h*i*(y_N - y_0)/PI + y_0);
        y_SW_err.push(h*i, abs(y_SW.v_c[i].second - y_A.v_c[i].second));
    }
    saveVectorPoint2DToFile(y_A.real(), "y_A.dat");
    saveVectorPoint2DToFile(y_SW_err.real(), "y_SW_err.dat");

    GnuplotPipe gp;
    //gp.sendLine(R"(plot "y_SW.dat" with lines, "y_A.dat" with lines)");
    gp.sendLine(R"(plot "y_SW_err.dat" with lines)");
}
