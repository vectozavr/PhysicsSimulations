//
// Created by Иван Ильин on 22.10.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

void crankNicolsonMethod(vector<vector<Point3D>>& map, double t_N = 5.0f, double h = 0.001, double tau  = 0.001) {
    map.clear();
    // the number of steps (calculated automatically)
    int N = 1 + 1.0f / h;
    int N_time = 1 + t_N / tau;

    // initial conditions
    for(int i = 0; i < N_time; i++) {
        map.push_back(vector<Point3D>(N));
        map[i][0].z = 0.0f;
        map[i][0].x = 0.0f;
        map[i][0].y = i*tau;

        map[i][N-1].z = 0.0f;
        map[i][N-1].x = 1.0f;
        map[i][N-1].y = i*tau;
    }
    for(int i = 0; i < N-1; i++) {
        //map[0][i].z = i * h * (1.0f - i * h) * (1.0f - i * h);
        map[0][i].z = 0.0f;
        map[0][i].y = 0.0f;
        map[0][i].x = i*h;
    }

    // for all time points m we should use sweep method:
    for(int m = 1; m < N_time; m++) {
        vector<vector<double>> matrix;
        vector<double> d(N-2);
        vector<double> y(N-2);
        // Av = d, here we fullfit the A matrix and d vector:
        for(int i = 0; i < N-2; i++) {
            matrix.push_back({});
            for(int j = 0; j < N-2; j++) {
                if(j == i-1) {
                    matrix[i].push_back(-0.5f*tau/(h*h));
                } else if(j == i) {
                    matrix[i].push_back(1.0f + tau/(h*h));
                } else if(j == i + 1) {
                    matrix[i].push_back(-0.5f*tau/(h*h));
                } else {
                    matrix[i].push_back(0.0f);
                }
            }
            d[i] = (map[m-1][i+1].z + 0.5f*tau*( 2.0f*i*h*(1.0f - i*h) + (map[m-1][i+2].z - 2.0f*map[m-1][i+1].z + map[m-1][i].z)/(h*h)));
        }

        // dv/dx(1, t) = 0:
        matrix[N-3][N-5] = 0.0f;
        matrix[N-3][N-4] = -1.0f;
        matrix[N-3][N-3] = 1.0f;
        d[N-3] = 0;

        // Newton' method
        for(int i = 1; i < N-2; i++) {
            double ksi = matrix[i][i-1] / matrix[i-1][i-1];
            matrix[i][i-1] = 0;
            matrix[i][i] -= ksi * matrix[i-1][i];
            d[i] -= ksi * d[i-1];
        }
        // back propagation
        y[N-3] = d[N-3] / matrix[N-3][N-3];
        for(int i = N-4; i >= 0; i--)
            y[i] = (d[i] - matrix[i][i+1] * y[i+1]) / matrix[i][i];

        for(int i = 0; i < y.size(); i++) {
            map[m][i+1].z = y[i];
            map[m][i+1].x = i*h;
            map[m][i+1].y = m*tau;
        }

        // dv/dx(1, t) = 0:
        map[m][y.size() + 1].z = y[y.size() - 1];

    }
}

int main() {

    // MOTION OF PLANETS
    ComplexPlot v_NM_err;
    ComplexPlot v_IM_err;
    ComplexPlot v_EX_err;

    // initial conditions
    double t_0 = 0.0f;
    double t_N = 0.6f;
    // integration steps
    double h    = 0.01;
    double tau  = 0.01;

    vector<vector<Point3D>> map;
    crankNicolsonMethod(map, t_N, h, tau);

    cout << "N_time = " << map.size() << endl;

    ComplexPlot maxV;
    for(int m = 0; m < map.size(); m++) {
        double d_maxV = map[m][0].z;
        for(int i = 1; i < map[m].size(); i++)
            if(d_maxV < map[m][i].z)
                d_maxV = map[m][i].z;
        maxV.push(tau*m, d_maxV);
    }

    saveVectorPoint3DToFile(map, "map.dat");
    saveVectorPoint2DToFile(maxV.real(), "maxV.dat");

    GnuplotPipe gp;
    //gp.sendLine(R"(plot "v_NM_err.dat" with lines, "v_IM_err.dat" with lines, "v_EX_err.dat" with lines)");

    gp.sendLine(R"(set xyplane relative 0)");
    gp.sendLine(R"(set view 60, 120)");

    gp.sendLine(R"(splot "map.dat" w pm3d)");
    //gp.sendLine(R"(plot "maxV.dat" with lines)");

}
