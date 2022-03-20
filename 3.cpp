#define _USE_MATH_DEFINES
#include <graphics.h>
#include <math.h>

#define SIZE 4

float coords[SIZE][2] = {{300, 309}, {304, 325}, {282, 311}, {276, 314}}; 

void line_b(int x1, int y1, int x2, int y2) {
    int dx = abs(x2 - x1), sx = x1 < x2 ? 1 : -1;
    int dy = abs(y2 - y1), sy = y1 < y2 ? 1 : -1; 
    int err = (dx > dy ? dx : -dy) / 2, e2;
    while (1) {
        putpixel(x1, y1, WHITE);
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

bool ccw(float a[2], float b[2], float c[2]) {
    return (c[1]-a[1])*(b[0]-a[0]) > (b[1]-a[1])*(c[0]-a[0]);
}

bool intersect(float a[2], float b[2], float c[2], float d[2]) {
    return (ccw(a, c, d) != ccw(b, c, d)) && (ccw(a, b, c) != ccw(a, b, d));
}

bool intersection(float line1[2][2], float line2[2][2], float *i_x, float *i_y) 
{
    float s1_x, s1_y, s2_x, s2_y;
    s1_x = line1[1][0] - line1[0][0];
    s1_y = line1[1][1] - line1[0][1];
    s2_x = line2[1][0] - line2[0][0];
    s2_y = line2[1][1] - line2[0][1];

    float s, t;
    s = (-s1_y * (line1[0][0] - line2[0][0]) + s1_x * (line1[0][1] - line2[0][1])) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (line1[0][1] - line2[0][1]) - s2_y * (line1[0][0] - line2[0][0])) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        if (i_x != NULL)
            *i_x = line1[0][0] + (t * s1_x);
        if (i_y != NULL)
            *i_y = line1[0][1] + (t * s1_y);
        return 1;
    }

    return 0;
}

bool in_poly(float x, float y, float coords[SIZE][2], float bounds[2][2]) {
    int intersections = 0;
    float pt_vec[2][2] = {{bounds[0][0]-1, y}, {x, y}};
    for (int i = 0; i < SIZE; i++)
    {
        float line[2][2] = {{coords[i][0]),             coords[i][1]},
                            {coords[(i + 1) % SIZE][0], coords[(i + 1) % SIZE][1]}};
        if (intersect(line[0], line[1], pt_vec[0], pt_vec[1]))
            intersections++;
    }
    return intersections & 1;
}

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

void fill(float coords[SIZE][2]) {
    float bounds[2][2] = {{coords[0][0], coords[0][1]}, {coords[0][0], coords[0][1]}};
    for (int i = 1; i < SIZE; i++) {
        int x = coords[i][0], y = coords[i][1];
        if (x < bounds[0][0])
            bounds[0][0] = x;
        else if (x > bounds[1][0])
            bounds[1][0] = x;
        if (y < bounds[0][1])
            bounds[0][1] = y;
        else if (y > bounds[1][1])
            bounds[1][1] = y;
    }
    int count = 0;
    for (int i = bounds[0][1]; i <= bounds[1][1]; i++) {
        bool found = false;
        int pts[20] = {0};
        for (int j = bounds[0][0]; j <= bounds[1][0]; j++) {
            if ((getpixel(j, i) && !getpixel(j + 1, i))) //|| (getpixel(j, i) && !getpixel(j - 1, i)))
            {
                // putpixel(j, i, BLUE);
                pts[count++] = j;
            }
            else if (getpixel(j, i) && !getpixel(j - 1, i))
            {
                // putpixel(j, i, GREEN);
                pts[count++] = j;
            }
        }
        if (count > 1)
        {
            for (int j = 0; j < count - 1; j++) {
                if (i == 310)
                    printf("%d ", in_poly(((pts[j] + pts[j+1]) / 2), i, coords, bounds));
                if (in_poly(((pts[j] + pts[j+1]) / 2), i, coords, bounds)) {
                    line(pts[j], i, pts[j + 1], i);    
                    // putpixel(((pts[j] + pts[j+1]) / 2), i, RED);
                }
                // putpixel(pts[j], i, BLUE);
            }
            if (i == 310)
                printf("\n");
        }
        count = 0;
    }
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
            line_b(coords[i][0], coords[i][1], coords[(i + 1) % SIZE][0], coords[(i + 1) % SIZE][1]);
        fill(coords);
        memset(offset, 0, sizeof(offset));
        angle = 0;
        swapbuffers();
    }
}