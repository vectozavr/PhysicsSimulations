//
// Created by vovam on 17.03.2020.
//

#include <SFML/Graphics.hpp>
#include <iostream>
#include <chrono>
#include "vemath.h"

using namespace std;
using namespace vemath;

[[nodiscard]] Point2D randomVelocity2D(double abs) {
    Point2D randomDir;

    double angle = 2*PI*(double)rand()/RAND_MAX;

    randomDir.x = abs*cos(angle);
    randomDir.y = abs*sin(angle);

    return randomDir;
}

double calculateShift(double elapsedTime, vector<Point2D>& turtles, double velocity, double R = 30*100, double dt = 0.01) {
    double time = 0;
    while (time < elapsedTime) {
        Point2D p = turtles[0] - turtles[1];
        if( p.abs() < 1 ) {
            return time;
        }
        for(int k = 0; k < turtles.size(); k++) {

            turtles[k].x += turtles[k].Vx*dt;
            turtles[k].y += turtles[k].Vy*dt;

            if(turtles[k].x > R || turtles[k].x < -R) {
                turtles[k].Vx *= -1;
                turtles[k].x += turtles[k].Vx*dt;
            }
            if(turtles[k].y > R || turtles[k].y < -R) {
                turtles[k].Vy *= -1;
                turtles[k].y += turtles[k].Vy*dt;
            }

            turtles[k].Vx += turtles[k].Ax * dt;
            turtles[k].Vy += turtles[k].Ay * dt;

            double abs = sqrt(turtles[k].Vx * turtles[k].Vx + turtles[k].Vy * turtles[k].Vy);
            if(abs > velocity) {
                turtles[k].Vx = velocity * turtles[k].Vx / abs;
                turtles[k].Vy = velocity * turtles[k].Vy / abs;
            }
        }
        time += dt;
    }
    return time;
}

double calculateTime(double R, double velocity, unsigned long long N) {
    return R/(velocity*cos(PI/2 - (double)PI/N));
}



int main() {
    srand(126);

    sf::RenderWindow window(sf::VideoMode(SCREEN_WIDTH, SCREEN_HEIGHT), "turtles");
    window.clear(sf::Color(255, 255, 255, 0));
    auto startTime = chrono::system_clock::now();
    auto tp1 = chrono::system_clock::now();
    auto tp2 = chrono::system_clock::now();

    unsigned long long N = 1000;
    double R = 300*100;
    double velocity = 100;
    double acceleration = 1.f;

    double r = 1000;

    int i_totalTime = 0;
    double totalTime = 0;

    int inflectedNumber = 0;

    vector<Point2D> turtles{N, {0,0}};
    for(int i = 0; i < turtles.size(); i++) {
        //turtles[i].x = R*cos((double)i*2*PI/N);
        //turtles[i].y = R*sin((double)i*2*PI/N);

        // In random positions:
        turtles[i].x = R*(-1 + 2*(double)rand()/RAND_MAX);
        turtles[i].y = R*(-1 + 2*(double)rand()/RAND_MAX);

        Point2D direction = randomVelocity2D(velocity);
        turtles[i].Vx = direction.x;
        turtles[i].Vy = direction.y;

        if(i == 0) {
            turtles[i].inflected = true;
            inflectedNumber++;
        }
    }

    turtles[0].x = 0;
    turtles[0].y = 0;

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

        double plus = calculateShift(1.f, turtles, velocity, R);
        totalTime += plus;

        for(int k = 0; k < turtles.size(); k++) {
            if(i_totalTime != (int)totalTime) {
                Point2D direction = randomVelocity2D(acceleration);

                turtles[k].Ax = velocity*direction.x;
                turtles[k].Ay = velocity*direction.y;

            }
        }

        //cout << totalTime << endl;

        //window.clear();     // отчистка

        window.clear(sf::Color(255, 255, 255, 0));
        sf::RectangleShape rectangle(sf::Vector2f(2*R*SCALE, 2*R*SCALE));
        rectangle.move(SCREEN_WIDTH/2 - R*SCALE, SCREEN_HEIGHT/2 - R*SCALE);
        rectangle.setFillColor(sf::Color(255, 255, 255));
        rectangle.setOutlineThickness(10.f); // устанавливаем толщину контура круга
        rectangle.setOutlineColor(sf::Color(80,220,50)); // устанавливаем цвет контура
        window.draw(rectangle);

        for(int i = 0; i < turtles.size(); i++) {

            int circleRadius = 2;

            sf::CircleShape circle(circleRadius);

            if(turtles[i].inflected) {
                circle.setFillColor(sf::Color(255, 0, 0));
                sf::CircleShape circle2(r*SCALE);
                circle2.setFillColor(sf::Color(255, 0, 0, 100));
                circle2.setPosition(SCREEN_WIDTH/2 + (int)turtles[i].x*SCALE - r*SCALE, SCREEN_HEIGHT/2 + (int)turtles[i].y*SCALE - r*SCALE);
                window.draw(circle2);

                if(i_totalTime != (int)totalTime) {
                    for(int k = 0; k < turtles.size(); k++) {
                        if(!turtles[k].inflected) {
                            if((turtles[i] - turtles[k]).abs() < r) {
                                if(((double)rand()/RAND_MAX) > 0.7) {
                                    if(!turtles[k].inflected)
                                        inflectedNumber++;
                                    turtles[k].inflected = true;

                                    //turtles[k].x = R*(-1 + 2*(double)rand()/RAND_MAX);
                                    //turtles[k].y = R*(-1 + 2*(double)rand()/RAND_MAX);
                                }
                            }
                        }
                    }
                }
            }
            else
                circle.setFillColor(sf::Color(0, 0, 0));

            circle.setPosition(SCREEN_WIDTH/2 + (int)turtles[i].x*SCALE - circleRadius/2, SCREEN_HEIGHT/2 + (int)turtles[i].y*SCALE - circleRadius/2);
            window.draw(circle);
        }
        i_totalTime = (int) totalTime;
        window.display();   // отображение

        cout << (int)totalTime << "  " << inflectedNumber << endl;

    }

    return 0;
}
