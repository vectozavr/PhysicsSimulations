//
// Created by ivan- on 09.05.2020.
//

#include <SFML/Graphics.hpp>
#include <iostream>
#include <chrono>
#include <cmath>

using namespace std;

struct Point2D {
    double x = 0;
    double y = 0;
};

int randomPoint(unsigned count) {
    return floor((double)count*rand()/RAND_MAX);
}

void nextPoint(double& x, double& y) {
    float nextX, nextY;
    /*
    float r = (float)rand()/RAND_MAX;
    if (r < 0.01) {         // Рост стебля
        nextX =  0;
        nextY =  0.16 * y;
    } else if (r < 0.86) {  // Рост маленьких листьев
        nextX =  0.85 * x + 0.04 * y;
        nextY = -0.0 * x + 0.85 * y + 1.6;
    } else if (r < 0.93) {  // Рост больших листьев справа
        nextX =  0.20 * x - 0.26 * y;
        nextY =  0.23 * x + 0.22 * y + 1.6;
    } else {                // Рост больших листьев слева
        nextX = -0.15 * x + 0.28 * y;
        nextY =  0.26 * x + 0.24 * y + 0.44;
    }
     */
    x = nextX;
    y = nextY;
}

int main() {
    cout << floor(2.9) << endl;

    int SCREEN_WIDTH    = 1080;
    int SCREEN_HEIGHT   = 1920;

    srand(12234);

    sf::ContextSettings settings;
    settings.antialiasingLevel = 4;

    sf::RenderWindow window(sf::VideoMode(SCREEN_WIDTH, SCREEN_HEIGHT), "point_frac", sf::Style::Default, settings);
    window.clear(sf::Color(255, 255, 255, 0));
    auto startTime = chrono::system_clock::now();
    auto tp1 = chrono::system_clock::now();
    auto tp2 = chrono::system_clock::now();

    double r0 = 3;
    unsigned long long N    = 3;
    double totalTime        = 0;
    double scale            = 100;

    double probability1 = 0.8;

    int r = 3;
    int R = 8;

    vector<Point2D> points;
    //points.push_back({0, 3});
    //points.push_back({-3.5, -1.5});
    //points.push_back({3, 2.8});
    //points.push_back({2, 1});
    //points.push_back({3, -2});
    //points.push_back({1.5, -3});
    for(int i = 0; i < N; i++) {
        points.push_back({r0*cos(i*2*3.1415/N), r0*sin(i*2*3.1415/N)});
    }

    Point2D cur = {(points[0].x + points[1].x)/2, (points[0].y + points[1].y)/2};

    for(int k = 0; k < points.size(); k++) {
        sf::CircleShape circle(R);
        circle.setFillColor(sf::Color(255, 0, 0, 155));
        circle.setPosition(SCREEN_WIDTH/2 + points[k].x*scale - R, SCREEN_HEIGHT/2 - points[k].y*scale - R);
        window.draw(circle);
    }


    int flag = 10000;

    while (window.isOpen())
    {
        tp2 = chrono::system_clock::now();
        chrono::duration <double> elapsedTime = tp2 - tp1;
        tp1 = tp2;
        double d_elapsedTime = elapsedTime.count();

        std::string title = "point_frac " + std::to_string((double)1/elapsedTime.count());
        window.setTitle(title);
        sf::Event event{};
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        //if(sf::Keyboard::isKeyPressed(sf::Keyboard::Space) && (flag > 15))
        {
            //window.clear(sf::Color(255, 255, 255, 0));

            flag = 0;
            sf::CircleShape circle(r);
            circle.setFillColor(sf::Color(255, 0, 0, 100));
            circle.setPosition(SCREEN_WIDTH / 2 + cur.x * scale - r,SCREEN_HEIGHT / 2 - cur.y * scale - r);
            window.draw(circle);

            int rand = randomPoint(points.size());
            //cout << points.size() << endl;
            Point2D vecTranslation = {(points[rand].x - cur.x) / 2, (points[rand].y - cur.y) / 2};

            cur.x += vecTranslation.x;
            cur.y += vecTranslation.y;
            //nextPoint(cur.x, cur.y);
            window.display();   // отображение
        }
        flag++;
        totalTime++;

    }

    system("pause");
    return 0;
}
