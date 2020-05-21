//
// Created by ivan- on 06.05.2020.
//
#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>

using namespace sf;
using namespace std;

const int WIDTH = 800, HEIGHT = 800, MAX_IT = 200, DEG = 1;
const double MAX_LEN = 2.0;
const int R = 255, G = 204, B = 204;
const int dR = 63, dG = 7, dB = 63;

double x1 = -1.6, __y1 = -1.6, x2 = 1.6, y2 = 1.6;
//double x1 = -0.1, __y1 = -0.1, x2 = 0.1, y2 = 0.1;

struct C_complex
{
    double a, b;
    C_complex(double _a = 0, double _b = 0)
    {
        a = _a;
        b = _b;
    }
    C_complex operator * (C_complex cur)
    {
        double c = cur.a;
        double d = cur.b;
        double _a = a * d + b * c;
        double _b = b * d - a * c;
        return C_complex(_a, _b);
    }
    C_complex operator + (C_complex cur)
    {
        return C_complex(a + cur.a, b + cur.b);
    }
    double get_len()
    {
        return sqrt(a * a + b * b);
    }
};

double get_x(int x)
{
    double res = (double) x * (x2 - x1);
    res /= WIDTH;
    return res + x1;
}

double get_y(int y)
{
    double res = (double) y * (y2 - __y1);
    res /= HEIGHT;
    return res + __y1;
}

void draw_point(sf::RenderWindow& window, int x, int y, int r, int g, int b){
    sf::Vertex line[] =
            {
                    sf::Vertex(sf::Vector2f(x, y), Color(r, g, b)),
                    sf::Vertex(sf::Vector2f(x, y + 1), Color(r, g, b))
            };
    window.draw(line, 2, sf::Lines);
}

void draw_fractal(sf::RenderWindow& window, C_complex c_0)
{
    for (int i = 0; i < WIDTH; ++i)
        for (int j = 0; j < HEIGHT; ++j)
        {
            int cnt = MAX_IT;
            double x = get_x(i);
            double y = get_y(j);
            C_complex z(y, x), c(y, x);
            for (int it = 0; it < MAX_IT; ++it)
            {
                z = z * z + c_0;
                if (z.get_len() >= MAX_LEN)
                {
                    cnt = it;
                    break;
                }
            }

            //int r = (dR - R) * cnt / MAX_IT + R;
            //int g = (dG - G) * cnt / MAX_IT + G;
            //int b = (dB - B) * cnt / MAX_IT + B;

            int r = 255 - 255*cnt / MAX_IT;
            int g = 255 - 255*cnt / MAX_IT;
            int b = 255 - 255*cnt / MAX_IT;

            draw_point(window, i, j, r, g, b);
        }
}

int main()
{
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "SFML works!", sf::Style::Default, settings);

    int iter = 38913000;

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();
        draw_fractal(window, {(double)iter/100000000*cos((double)iter/1000000), (double)iter/100000000*sin((double)iter/1000000)});
        window.display();
        Vector2i position = Mouse::getPosition(window);
        Vector2i position1;
        //while (Mouse::isButtonPressed(Mouse::Left))
        //{
        //    position1 = Mouse::getPosition(window);
        //}
        //double _x1 = get_x(position.x);
        //double _y1 = get_y(position.y);
        //double _x2 = get_x(position1.x);
        //double _y2 = get_y(position1.y);
        //x1 = min(_x1, _x2);
        //__y1 = min(_y1, _y2);
        //x2 = max(_x1, _x2);
        //y2 = max(_y1, _y2);
        iter++;

    }
    cout << iter << endl;
    return 0;
}
