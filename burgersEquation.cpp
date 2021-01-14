//
// Created by Иван Ильин on 13.01.2021.
//

#include <cmath>
#include <SFML/Graphics.hpp>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

void explicitMethod(vector<vector<Point3D>>& map, double e, double t_N = 5.0f, double h = 0.001, double tau  = 0.001) {
    map.clear();
    // the number of steps (calculated automatically)
    int N = 1 + 10.0f / h;
    int N_time = 1 + t_N / tau;

    // initial conditions
    for(int i = 0; i < N_time; i++) {
        map.push_back(vector<Point3D>(N));
        map[i][0].z = 0;
        map[i][0].x = 0;
        map[i][0].y = i*tau;

        map[i][N-1].z = 0;
        map[i][N-1].x = 10;
        map[i][N-1].y = i*tau;
    }
    for(int i = 0; i < N-1; i++) {
        if ((i*h >= 1) && (i*h <= 2))
            map[0][i].z = 1.0;
        else if ((i*h >= 3) && (i*h <= 5))
            map[0][i].z = 2.0;
        else
            map[0][i].z = 0.0;

        map[0][i].y = 0;
        map[0][i].x = i*h;
    }

    for(int m = 1; m < N_time; m++) {
        for(int j = 1; j < N-1; j++) {
            map[m][j].x = j*h;
            map[m][j].y = m*tau;
            map[m][j].z = tau*e*(map[m-1][j+1].z - 2*map[m-1][j].z + map[m-1][j-1].z)/(h*h) - tau*map[m-1][j].z*(map[m-1][j+1].z - map[m-1][j-1].z)/(2.0*h) + map[m-1][j].z;
        }
    }

}

void crankNicolsonMethod(vector<vector<Point3D>>& map, double e, double t_N = 10.0f, double h = 0.001, double tau  = 0.001) {
    map.clear();
    // the number of steps (calculated automatically)
    int N = 1 + 10.0f / h;
    int N_time = 1 + t_N / tau;

    // initial conditions
    for(int i = 0; i < N_time; i++) {
        map.push_back(vector<Point3D>(N));
        map[i][0].z = 0;
        map[i][0].x = 0;
        map[i][0].y = i*tau;

        map[i][N-1].z = 0;
        map[i][N-1].x = 10;
        map[i][N-1].y = i*tau;
    }
    for(int i = 0; i < N-1; i++) {
        if ((i*h >= 1) && (i*h <= 2))
            map[0][i].z = 1.0;
        else if ((i*h >= 3) && (i*h <= 5))
            map[0][i].z = 2.0;
        else
            map[0][i].z = 0.0;

        map[0][i].y = 0;
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
                    matrix[i].push_back(-map[m-1][j].z/(4.0*h) - e/(2.0*h*h));
                } else if(j == i) {
                    matrix[i].push_back(1.0f/tau + e/(h*h));
                } else if(j == i + 1) {
                    matrix[i].push_back(map[m-1][j].z/(4.0*h) - e/(2.0*h*h));
                } else {
                    matrix[i].push_back(0.0f);
                }
            }
            d[i] = (map[m-1][i+1].z/tau - map[m-1][i+1].z/(4.0*h)*(map[m-1][i+2].z - map[m-1][i].z) + e/(2.0*h*h)*(map[m-1][i+2].z - 2.0*map[m-1][i+1].z + map[m-1][i].z));
        }

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
    }
}

int main() {
    // initial conditions
    double t_0 = 0.0f;
    double t_N = 10.0f;
    // integration steps
    double h    = 0.1;
    double tau  = 0.1;

    double e = 0.0;
    double de = 0.1;

    vector<vector<Point3D>> map_e;
    vector<vector<Point3D>> map;

    while (e <= 1.0) {
        crankNicolsonMethod(map, e, t_N, h, tau);

        map_e.push_back(map.back());
        for(auto& elem : map_e.back())
            elem.y = e;
        e += de;
    }

    crankNicolsonMethod(map, 0.5, t_N, 0.05, 0.05);
    saveVectorPoint3DToFile(map, "map.dat");
    saveVectorPoint3DToFile(map_e, "map_e.dat");

    GnuplotPipe gp;

    gp.sendLine(R"(set xyplane relative 0)");
    gp.sendLine(R"(set view 60, 120)");

    gp.sendLine(R"(splot "map_e.dat" w pm3d)");
}

