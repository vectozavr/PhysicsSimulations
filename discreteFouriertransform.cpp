//
// Created by ivan- on 11.02.2020.
//
#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>

#include "settings.h"
#include "gnuplot.h"

struct Complex {
    double real = 0;
    double imagine = 0;
};

struct Point2D {
    double x = 0;
    double y = 0;

    Point2D& operator+=(const Point2D& point2D) { this->x += point2D.x; this->y += point2D.y; }
    Point2D& operator=(const Point2D& point2D) { this->x = point2D.x; this->y = point2D.y; return *this; }
    Point2D& operator*(double number) { this->x *= number; this->y *= number; }
    Point2D operator-(const Point2D& point2D) const { return {this->x - point2D.x, this->y - point2D.y}; }
    Point2D operator+(const Point2D& point2D) const { return {this->x + point2D.x, this->y + point2D.y}; }

    Point2D normalize() { return Point2D{this->x/abs(), this->y/abs()};}
    double abs() {return sqrt(x*x + y*y); }
};

bool saveVectorPoint2DToFile(const std::vector<Point2D>& data, const std::string& fileName = "data.dat") {
    std::ofstream _ofstream(fileName);
    if(!_ofstream.is_open())
        return false;

    for(auto d : data)
        _ofstream << d.x << "   " << d.y << std::endl;

    _ofstream.close();
    return true;
}

void fourierTransform(const std::vector<Point2D>& data, std::vector<Point2D>& transform) {
    transform.clear();
    for(int k = 0; k < data.size(); k++) {
        Point2D transformed = {0, 0};
        for(int n = 0; n < data.size(); n++) {
            double p = 2 * PI * k * n / data.size();

            transformed.x += data[n].y * cos(p);
            transformed.y -= data[n].y * sin(p);
        }
        transform.push_back({2 * PI * k / data.size(), transformed.abs() / data.size()});
    }
}

void addNoice(std::vector<Point2D>& data, double noiseAmplitude = 0.01) {
    srand(time(NULL));
    for(auto d : data) {
        //d.x += rand()/RAND_MAX;
        //d.y +=
    }
}

int main() {
    std::vector<Point2D> data;
    std::vector<Point2D> transform;

    int quantity = 1000;
    for(int i = -quantity; i < quantity; i++) {
        if(i < -quantity / 10 || i > quantity / 10)
            data.push_back({(double)i/quantity*10, 0});
        else
            data.push_back({(double)i/quantity*10, 1});
    }
    //for(int i = -quantity; i < quantity; i++) {
    //    double x = (double)i/quantity * 2 * PI;
    //    data.push_back({x, sin(x) + sin(5*x)/2});
    //}

    fourierTransform(data, transform);



    saveVectorPoint2DToFile(data, "data.dat");
    saveVectorPoint2DToFile(transform, "transform.dat");

    GnuplotPipe gp;
    gp.sendLine(R"(set xrange [-10:10])");
    gp.sendLine(R"(set yrange [-2:2])");
    gp.sendLine(R"(plot "data.dat" with lines, "transform.dat" with lines)");

    return 0;

}