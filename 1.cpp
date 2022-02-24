#define _USE_MATH_DEFINES
#include <graphics.h>
#include <math.h>

float coords[2][4] = {{150, 150, 450, 150}, {150, 250, 450, 250}}; 

void translate(float *coords, float *offset) {
    coords[0] += offset[0];
    coords[1] += offset[1];
    coords[2] += offset[0];
    coords[3] += offset[1];
}

void rotate(float *coords, int angle) {
    double theta = (double)(angle % 180) * M_PI / 180;
    float x = coords[2];
    float y = coords[3];
    x = coords[2] * cos(theta) + coords[3] * sin(theta); 
    y = -coords[2] * sin(theta) + coords[3] * cos(theta); 
    coords[2] = x;
    coords[3] = y;
}

void scale(float *coords, float prev_coef, float coef) {
    for (int i = 0; i < 4; i++)
        coords[i] = coords[i] / prev_coef * coef;
}

void line_b(int x1, int y1, int x2, int y2) {
    int dx = abs(x2 - x1), sx = x1 < x2 ? 1 : -1;
    int dy = abs(y2 - y1), sy = y1 < y2 ? 1 : -1; 
    int err = (dx > dy ? dx : -dy) / 2, e2;
    while (1) {
        putpixel(x1, y1, CYAN);
        if (x1 == x2 && y1 == y2)
            break;
        e2 = err;
        if (e2 >-dx) {
            err -= dy;
            x1 += sx;
        }
        if (e2 < dy) {
            err += dx;
            y1 += sy;
        }
    }
}
  
int main()
{
    float offset[2] = {0, 0};
    int angle = 0;
    int speed = 5;
    float coef = 1, prev_coef = 1;
    float *cur_coords = coords[0];
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
            else if (ch == '1') 
                cur_coords = coords[0];
            else if (ch == '2') 
                cur_coords = coords[1];
            else if (ch == '+')
                speed += 1;
            else if (ch == '-')
                if (speed > 1)
                    speed -= 1;
        }
        cleardevice();
        translate(cur_coords, offset);
        float origin[2] = {(cur_coords[0]+cur_coords[2])/2, (cur_coords[1]+cur_coords[3])/2};
        offset[0] = -cur_coords[0];
        offset[1] = -cur_coords[1];
        translate(cur_coords, offset);
        if (prev_coef != coef) {
            scale(cur_coords, prev_coef, coef);
            prev_coef = coef;
        }
        rotate(cur_coords, angle);
        float center[2] = {-(cur_coords[0]+cur_coords[2])/2, -(cur_coords[1]+cur_coords[3])/2};
        translate(cur_coords, center);
        translate(cur_coords, origin);
        line(coords[0][0], coords[0][1], coords[0][2], coords[0][3]);
        line_b(coords[1][0], coords[1][1], coords[1][2], coords[1][3]);
        memset(offset, 0, sizeof(offset));
        angle = 0;
        swapbuffers();
    }
}