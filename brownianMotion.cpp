//
// Created by ivan- on 03.05.2020.
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


int main() {
    srand(126);

    sf::ContextSettings settings;
    settings.antialiasingLevel = 4;

    sf::RenderWindow window(sf::VideoMode(SCREEN_WIDTH, SCREEN_HEIGHT), "turtles", sf::Style::Default, settings);
    window.clear(sf::Color(255, 255, 255, 0));
    auto startTime = chrono::system_clock::now();
    auto tp1 = chrono::system_clock::now();
    auto tp2 = chrono::system_clock::now();

    unsigned long long N = 1;
    double R = 500*100;
    double velocity = 500;
    double acceleration = 1.f;

    int i_totalTime = 0;
    double totalTime = 0;

    int steps = 0;

    double scale = 0.01;

    vector<Point2D> turtles{N, {0,0}};
    for(int i = 0; i < turtles.size(); i++) {

        // In random positions:
        turtles[i].x = R*(-1 + 2*(double)rand()/RAND_MAX);
        turtles[i].y = R*(-1 + 2*(double)rand()/RAND_MAX);

        Point2D direction = randomVelocity2D(velocity);
        turtles[i].Vx = direction.x;
        turtles[i].Vy = direction.y;
    }

    turtles[0].x = -2500;
    turtles[0].y =  2500;

    vector<Point2D> history;

    for(int k = 0; k < 100000; k++) {
        double plus = calculateShift(0.05f, turtles, velocity, R);
        totalTime += plus;

        for(int k = 0; k < turtles.size(); k++) {
            if(i_totalTime != (int)totalTime) {
                Point2D direction = randomVelocity2D(acceleration);

                turtles[k].Vx = velocity*direction.x;
                turtles[k].Vy = velocity*direction.y;
            }
        }

        history.push_back({turtles[0].x, turtles[0].y});
    }

    while (window.isOpen())
    {
        vector<Point2D> lastPositions = turtles;

        if(scale < 1)
            scale *= 1.01;

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

        //cout << totalTime << endl;
        //window.clear();     // отчистка

        window.clear(sf::Color(255, 255, 255, 0));
        sf::RectangleShape rectangle(sf::Vector2f(2*R*SCALE, 2*R*SCALE));
        rectangle.move(SCREEN_WIDTH/2 - R*SCALE, SCREEN_HEIGHT/2 - R*SCALE);
        rectangle.setFillColor(sf::Color(255, 255, 255));
        rectangle.setOutlineThickness(10.f); // устанавливаем толщину контура круга
        rectangle.setOutlineColor(sf::Color(80,220,50)); // устанавливаем цвет контура
        //window.draw(rectangle);

        for(int i = 0; i < turtles.size(); i++) {

            int circleRadius = 1;

            sf::CircleShape circle(circleRadius);

            circle.setFillColor(sf::Color(255, 255, 255, 50));

            circle.setPosition(SCREEN_WIDTH/2 + (int)turtles[i].x*SCALE - circleRadius/2, SCREEN_HEIGHT/2 + (int)turtles[i].y*SCALE - circleRadius/2);
            //window.draw(circle);

            for(int s = steps; s < history.size()-1; s++) {
                sf::Vertex line[] =
                        {
                                sf::Vertex(sf::Vector2f(SCREEN_WIDTH / 2 + (int) history[s].x * scale,
                                                        SCREEN_HEIGHT / 2 + (int) history[s].y * scale)),
                                sf::Vertex(sf::Vector2f(SCREEN_WIDTH / 2 + (int) history[s+1].x * scale,
                                                        SCREEN_HEIGHT / 2 + (int) history[s+1].y * scale))
                        };

                line->color = sf::Color{0, 0, 0, 50};

                window.draw(line, 2, sf::Lines);
            }
            //steps = history.size()-1;

        }
        i_totalTime = (int) totalTime;
        window.display();   // отображение

    }

    system("pause");

    return 0;
}
