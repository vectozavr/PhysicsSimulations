//
// Created by ivan- on 17.02.2020.
//

#ifndef PHYSICSSIMULATIONS_VEMATH_H
#define PHYSICSSIMULATIONS_VEMATH_H

#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include "settings.h"

namespace vemath {
    // struct Point2D allows you to do simple calculations with vectors really fast and comfortably.
    struct Point2D {
        double x = 0;
        double y = 0;

        Point2D &operator+=(const Point2D &point2D) {
            this->x += point2D.x;
            this->y += point2D.y;
        }
        Point2D &operator=(const Point2D &point2D) {
            this->x = point2D.x;
            this->y = point2D.y;
            return *this;
        }
        Point2D &operator*(double number) {
            this->x *= number;
            this->y *= number;
        }
        double operator*(const Point2D &point2D) const { return this->x * point2D.x + this->y * point2D.y; }
        Point2D operator-(const Point2D &point2D) const { return {this->x - point2D.x, this->y - point2D.y}; }
        Point2D operator+(const Point2D &point2D) const { return {this->x + point2D.x, this->y + point2D.y}; }
        Point2D normalize() { return Point2D{this->x / abs(), this->y / abs()}; }
        double abs() { return sqrt(x * x + y * y); }
    };
    // struct Complex allows you to store complex numbers.
    struct Complex {
        double real = 0;
        double imagine = 0;

        Complex operator*(const Complex &complex) const { return {real*complex.real - imagine*complex.imagine, imagine*complex.real + real*complex.imagine}; }
        Complex &operator+=(const Complex &complex) {
            this->real += complex.real;
            this->imagine += complex.imagine;
        }

        [[nodiscard]] double abs() const {return sqrt(real*real + imagine*imagine); }
    };
    struct ComplexPlot {
        std::vector<std::pair<double, Complex>> v_c;

        [[nodiscard]] std::vector<Point2D> abs() const{
            std::vector<Point2D> result;
            for(int i = 0; i < v_c.size(); i++)
                result.push_back({v_c[i].first, v_c[i].second.abs()});
            return result;
        }

        void push(double x, Complex z){ v_c.emplace_back(x, z); }
    };

    // save 2D plot from <data> to "filename"
    bool saveVectorPoint2DToFile(const std::vector<Point2D> &data, const std::string &fileName = "data.dat");
    // simple fourier transform of 2D plot. Without imagine component.
    void fourierTransform(const ComplexPlot& data, ComplexPlot& transform);
    // add some noise to 2D plot <data> with amplitude <noiseAmplitude>
    void inverseFourierTransform(const ComplexPlot& data, ComplexPlot& transform);
    // add some noise to 2D plot <data> with amplitude <noiseAmplitude>
    void addNoise(ComplexPlot& data, double noiseAmplitude = 0.1);

}
#endif //PHYSICSSIMULATIONS_VEMATH_H
