//
// Created by Иван Ильин on 15.10.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

int main() {

    // MOTION OF PLANETS
    ComplexPlot psi;
    ComplexPlot psi_A;
    ComplexPlot psi_err;

    // range of integration
    double x_0 = -5;
    double x_N = 5;
    // integration step
    double h = 0.01;

    double deltaMin = 0.001f;
    double delta = 1.0f;

    // the number of steps (calculated automatically)
    int N = (x_N - x_0) / h;

    // initial lambda = 1 and psi = (1, 1, .., 1)
    double lambda = 0.0f;
    vector<double> v_psi(N, 1.0f/N);

    while(delta > deltaMin) {
        vector<vector<double>> matrix;
        vector<double> d(N);
        vector<double> y(N);
        // (A - lE)y = d, here we fullfit the A matrix and d vector
        for(int i = 0; i < N; i++) {
            matrix.push_back({});
            for(int j = 0; j < N; j++) {
                if(j == i-1) {
                    matrix[i].push_back(-1.0f/(2.0f*h*h));
                } else if(j == i) {
                    matrix[i].push_back(1.0f/(h*h) + (x_0 + i*h)*(x_0 + i*h)/2.0f);
                } else if(j == i + 1) {
                    matrix[i].push_back(-1.0f/(2.0f*h*h));
                } else {
                    matrix[i].push_back(0.0f);
                }
            }
            d[i] = (v_psi[i]);
        }

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

        double norm2 = 0;
        double norm1 = 0;

        for(int i = 0; i < N; i++) {
            norm1 += v_psi[i] * v_psi[i];
            norm2 += y[i] * y[i];

            v_psi[i] = y[i];
        }

        delta = abs(lambda - sqrt(norm1) / sqrt(norm2));
        lambda = sqrt(norm1) / sqrt(norm2);
    }
    double nor = 0;
    double max = v_psi[0];
    for(int i = 0; i < N; i++) {
        nor += v_psi[i] * v_psi[i];
        if(max < v_psi[i])
            max = v_psi[i];
    }
    for(int i = 0; i < N; i++) {
        v_psi[i] = v_psi[i] / max;
    }

    // analytical solution
    for(int i = 0; i < N; i++) {
        psi.push(x_0 + h*i, v_psi[i]);
        psi_A.push(x_0 + h*i, exp(-(x_0+i*h)*(x_0+i*h)/2.0f));
        psi_err.push(x_0 + h*i, abs(psi_A.v_c[i].second - v_psi[i]));
    }

    saveVectorPoint2DToFile(psi.real(), "psi.dat");
    saveVectorPoint2DToFile(psi_A.real(), "psi_A.dat");
    saveVectorPoint2DToFile(psi_err.real(), "psi_err.dat");

    cout << lambda << endl;

    GnuplotPipe gp;
    //gp.sendLine(R"(plot "psi_err.dat" with lines)");
    gp.sendLine(R"(plot "psi_A.dat" with lines, "psi.dat" with lines )");
}
