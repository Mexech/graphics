#define _USE_MATH_DEFINES
#include <graphics.h>
#include <math.h>

#define SIZE 4

float coords[SIZE][2] = {{150, 150}, {450, 150}, {450, 300}, {150, 250}}; 

void translate(float *coord, float *offset) {
    coord[0] += offset[0];
    coord[1] += offset[1];
}

void rotate(float *coord, int angle) {
    double theta = (double)(angle % 180) * M_PI / 180;
    float x = coord[0], y = coord[1];
    x = coord[0] * cos(theta) + coord[1] * sin(theta); 
    y = -coord[0] * sin(theta) + coord[1] * cos(theta); 
    coord[0] = x, coord[1] = y;
}

void scale(float *coord, float prev_coef, float coef) {
    for (int i = 0; i < 2; i++)
        coord[i] = coord[i] / prev_coef * coef;
}
  
int main()
{
    float offset[2] = {0, 0};
    int angle = 0;
    int speed = 5;
    float coef = 1, prev_coef = 1;
    int gd = DETECT, gm;
    initgraph(&gd, &gm, (char *)"");
    while (1) {
        if (kbhit()) {
            char ch = getch();
            if (ch == '\n')
                exit(0);
            else if (ch == 'w')
                offset[1] -= speed;
            else if (ch == 's')
                offset[1] += speed;
            else if (ch == 'a')
                offset[0] -= speed;
            else if (ch == 'd')
                offset[0] += speed;
            else if (ch == 'e')
                angle += speed;
            else if (ch == 'q')
                angle -= speed;
            else if (ch == 'x') {
                prev_coef = coef;
                coef += 0.01 * speed;
            }
            else if (ch == 'z') {
                prev_coef = coef;
                if (coef > 0.01 * speed)
                    coef -= 0.01 * speed;
            }
            else if (ch == '+')
                speed += 1;
            else if (ch == '-')
                if (speed > 1)
                    speed -= 1;
        }
        cleardevice();
        
        float origin[2] = {0};
        for (int i = 0; i < SIZE; i++)
            for (int j = 0; j < 2; j++)
                origin[j] += coords[i][j];
        origin[0] /= -SIZE; origin[1] /= -SIZE;
        for (int i = 0; i < SIZE; i++) {
            translate(coords[i], offset);
            translate(coords[i], origin);
            rotate(coords[i], angle);
            scale(coords[i], prev_coef, coef);
        }
        prev_coef = coef;
        origin[0] *= -1; origin[1] *= -1;
        for (int i = 0; i < SIZE; i++)
            translate(coords[i], origin);

        for (int i = 0; i < SIZE; i++)
            line(coords[i][0], coords[i][1], coords[(i + 1) % SIZE][0], coords[(i + 1) % SIZE][1]);
        memset(offset, 0, sizeof(offset));
        angle = 0;
        swapbuffers();
    }
}