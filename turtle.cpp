#include <SFML/Graphics.hpp>
#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>

#include "settings.h"

using namespace std;

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

double calculateShift(double elapsedTime, vector<Point2D>& turtles, double velocity, double dt = 0.001) {
    double time = 0;
    while (time < elapsedTime) {
        Point2D p = turtles[0] - turtles[1];
        if( p.abs() < 1 ) {
            return time;
        }
        for(int k = 0; k < turtles.size(); k++) {
            int m = k != turtles.size() - 1 ? k+1 : 0;
            Point2D direction = (turtles[m] - turtles[k]).normalize();

            turtles[k].x += velocity*dt*direction.x;
            turtles[k].y += velocity*dt*direction.y;
        }
        time += dt;
    }
    return time;
}

double calculateTime(double R, double velocity, unsigned long long N) {
    return R/(velocity*cos(PI/2 - (double)PI/N));
}



int main() {
    sf::RenderWindow window(sf::VideoMode(SCREEN_WIDTH, SCREEN_HEIGHT), "turtles");

    auto startTime = chrono::system_clock::now();
    auto tp1 = chrono::system_clock::now();
    auto tp2 = chrono::system_clock::now();

    unsigned long long N = 50;
    double R = 30*100;
    double velocity = 3;

    double totalTime = 0;

    vector<Point2D> turtles{N, {0,0}};
    for(int i = 0; i < turtles.size(); i++) {
        //turtles[i].x = R*cos((double)i*2*PI/N);
        //turtles[i].y = R*sin((double)i*2*PI/N);

        // In random positions:
        turtles[i].x = R*(-1 + 2*(double)rand()/RAND_MAX);
        turtles[i].y = R*(-1 + 2*(double)rand()/RAND_MAX);
    }

    double thTime = calculateTime(R, velocity, N);
    //double simTime = calculateShift(10000, turtles, velocity);

    while (window.isOpen())
    {
        vector<Point2D> lastPositions = turtles;


        tp2 = chrono::system_clock::now();
        chrono::duration <double> elapsedTime = tp2 - tp1;
        tp1 = tp2;
        double d_elapsedTime = elapsedTime.count();

        std::string title = "turtles " + std::to_string((double)1/elapsedTime.count());
        window.setTitle(title);
        sf::Event event{};
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        double plus = calculateShift(1.0f, turtles, velocity);
        totalTime += plus;

        //cout << totalTime << endl;

        //window.clear();     // отчистка
        for(int i = 0; i < turtles.size(); i++) {
            sf::CircleShape circle(1.f);
            circle.setFillColor(sf::Color(255, 125, 125));
            circle.setPosition(SCREEN_WIDTH/2 + (int)turtles[i].x*SCALE, SCREEN_HEIGHT/2 + (int)turtles[i].y*SCALE);
            window.draw(circle);

        }
        window.display();   // отображение
    }

    cout << "TheoTime = " << thTime << endl;
    cout << "CalcTime = " << totalTime << endl;
    cout << "Precision with " << 100*( 1 - abs(thTime - totalTime)/thTime ) << "%" << endl;

    return 0;
}