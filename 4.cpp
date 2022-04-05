#define _USE_MATH_DEFINES
#include <graphics.h>
#include <math.h>

#define COORDS_SIZE 5
#define GRAPH_SIZE 8
#define INIT_SPEED 1

float deg2rad = M_PI / 180.0;

float coords[COORDS_SIZE][3] = {{0.20, 0.5, 3.0},
                                {0.00, 2, 2.0},
                                {1.00, 2, 2.0},
                                {1.00, 2, 3.0},
                                {-1.0, 2, 3.0}};

int graph[GRAPH_SIZE][2] = {{1, 0}, {1, 4}, {1, 2}, {3, 0}, {3, 2}, {3, 4}, {2, 0}, {4, 0}};

void translate(float *coord, float *offset) {
    for (int i = 0; i < 3; i++)
        coord[i] += offset[i];
}

void rotate(float *coord, float *angles) {
    float x = coord[0], y = coord[1], z = coord[2];
    float c[3] = {cos(angles[0]*deg2rad), cos(angles[1]*deg2rad), cos(angles[2]*deg2rad)};
    float s[3] = {sin(angles[0]*deg2rad), sin(angles[1]*deg2rad), sin(angles[2]*deg2rad)};
    coord[0] = c[1]*(s[2]*y + c[2]*x) - s[1]*z;
    coord[1] = s[0]*(c[1]*z + s[1]*(s[2]*y + c[2]*x)) + c[0]*(c[2]*y - s[2]*x);
    coord[2] = c[0]*(c[1]*z + s[1]*(s[2]*y + c[2]*x)) - s[0]*(c[2]*y - s[2]*x);
}

void scale(float *coord, float prev_coef, float coef) {
    for (int i = 0; i < 3; i++)
        coord[i] = coord[i] / prev_coef * coef;
}

void drawMesh(float coords[][3], int graph[][2]) {
    for (int i = 0; i < GRAPH_SIZE; i++) {
        float *pt1 = coords[graph[i][0]];
        float *pt2 = coords[graph[i][1]];
        line(pt1[0], pt1[1], pt2[0], pt2[1]);
    }
}

void mult(float* pt, float* mat, float* ret) {
    for (int i = 0; i < 16; i++)
        ret[i / 4] += (i % 4 != 3 ? pt[i % 4] : 1)*mat[i];
}

void PerspectiveFOV(float fov, float aspect, float near_c, float far_c, float* ret) {
    float yScale = 1/tan(deg2rad*fov/2);
    float xScale = yScale/aspect;
    float diff = near_c - far_c;
    float m[] = {
        xScale, 0, 0, 0,
        0, yScale, 0, 0,
        0, 0, (far_c + near_c)/diff, (2*near_c*far_c)/diff, 
        0, 0, 1, 0
    };
    memcpy(ret, m, sizeof(float)*16);
}

int main() {
    float offset[3] = {0, 0, 0}, angles[3] = {0, 0, 0};
    int angle = 0;
    float speed = INIT_SPEED;
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
                offset[2] += speed;
            else if (ch == 'q') 
                offset[2] -= speed;
            else if (ch == '6') 
                angles[2] -= speed;
            else if (ch == '4') 
                angles[2] += speed;
            else if (ch == '8') 
                angles[0] -= speed;
            else if (ch == '2') 
                angles[0] += speed;
            else if (ch == '9') 
                angles[1] -= speed;
            else if (ch == '7') 
                angles[1] += speed;
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

        float origin[3] = {0}, neg_origin[3] = {0};
        for (int i = 0; i < COORDS_SIZE; i++)
            for (int j = 0; j < 3; j++)
                origin[j] += coords[i][j];
        for (int i = 0; i < 3; i++) {
            origin[i] /= COORDS_SIZE;
            neg_origin[i] = -origin[i];
        }

        float projection_mat[16] = {0};
        PerspectiveFOV(90, getwindowwidth()/getwindowheight(), 0.1, 100, projection_mat);
        float projected_coords[COORDS_SIZE][3];
        for (int i = 0; i < COORDS_SIZE; i++) {
            float projected_pt[4] = {0};
            translate(coords[i], offset);
            translate(coords[i], neg_origin);
            rotate(coords[i], angles);
            scale(coords[i], prev_coef, coef);
            translate(coords[i], origin);
            mult(coords[i], projection_mat, projected_pt);
            for (int j = 0; j < 3; j++)
                projected_coords[i][j] = 0.5 * (projected_pt[j] / projected_pt[3] + 1);
            projected_coords[i][0] *= getwindowwidth();
            projected_coords[i][1] *= getwindowheight();
        }

        drawMesh(projected_coords, graph);
        swapbuffers();

        prev_coef = coef;
        memset(offset, 0, sizeof(offset));
        memset(angles, 0, sizeof(angles));
    }
    return 0;
}