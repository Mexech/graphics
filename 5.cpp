#define _USE_MATH_DEFINES
#include <graphics.h>
#include <math.h>

#define COORDS_SIZE 5
#define INIT_SPEED 1
#define TRIS_SIZE 18
#define FACES_SIZE TRIS_SIZE / 3

float deg2rad = M_PI / 180.0;

float coords[COORDS_SIZE][3] = {{0.20, 0.5, 3.0},
                                {0.00, 2, 2.0},
                                {0.50, 2, 2.0},
                                {1.00, 2, 3.0},
                                {-1.0, 2, 3.0}};

int tris[TRIS_SIZE] = {4, 0, 1,
                       1, 0, 2,
                       2, 0, 3,
                       3, 0, 4,
                       1, 2, 4,
                       4, 2, 3};

int clrs[FACES_SIZE] = {BLUE, GREEN, RED, MAGENTA, CYAN, CYAN};

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

void drawMesh(float coords[][3], int tris[]) {
    for (int i = 0; i < TRIS_SIZE; i += 3) {
        for (int j = 0; j < 3; j++) {
            float *pt1 = coords[tris[i + j]];
            float *pt2 = coords[tris[i + (j + 1) % 3]];
            line(pt1[0], pt1[1], pt2[0], pt2[1]);
        }
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

int cmp(const void * a, const void * b) {
   return ( *(float*)a - *(float*)b );
}

bool ccw(float a[2], float b[2], float c[2]) {
    return (c[1]-a[1])*(b[0]-a[0]) > (b[1]-a[1])*(c[0]-a[0]);
}

bool intersect(float a[2], float b[2], float c[2], float d[2]) {
    return (ccw(a, c, d) != ccw(b, c, d)) && (ccw(a, b, c) != ccw(a, b, d));
}

float dot(float *a, float *b) {
    float res = 0;
    for (int i = 0; i < 3; i++)
        res += a[i]*b[i];
    return res;
}

float cross(float *a, float *b, float *c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

bool intersection(float line1[2][2], float line2[2][2], float *i_x, float *i_y) 
{
    float c1_x, c1_y, c2_x, c2_y;
    c1_x = line1[1][0] - line1[0][0];
    c1_y = line1[1][1] - line1[0][1];
    c2_x = line2[1][0] - line2[0][0];
    c2_y = line2[1][1] - line2[0][1];

    float s, t;
    s = (-c1_y * (line1[0][0] - line2[0][0]) + c1_x * (line1[0][1] - line2[0][1])) / (-c2_x * c1_y + c1_x * c2_y);
    t = ( c2_x * (line1[0][1] - line2[0][1]) - c2_y * (line1[0][0] - line2[0][0])) / (-c2_x * c1_y + c1_x * c2_y);
    
    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        if (i_x != NULL)
            *i_x = line1[0][0] + (t * c1_x);
        if (i_y != NULL)
            *i_y = line1[0][1] + (t * c1_y);
        return 1;
    }

    return 0;
}

bool in_poly(float x, float y, float coords[3][2], float bounds[2][2]) {
    int intersections = 0;
    float pt_vec[2][2] = {{bounds[0][0] - 1, y}, {x, y}};
    for (int i = 0; i < 3; i++)
    {
        float line[2][2] = {{coords[i][0],           coords[i][1]},
                            {coords[(i + 1) % 3][0], coords[(i + 1) % 3][1]}};
        if (intersect(line[0], line[1], pt_vec[0], pt_vec[1]))
            intersections++;
    }
    return intersections & 1;
}

void fill(float coords[3][2]) {
    float bounds[2][2] = {{coords[0][0], coords[0][1]}, {coords[0][0], coords[0][1]}};
    for (int i = 1; i < 3; i++) {
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
    for (int i = bounds[0][1]; i <= bounds[1][1]; i++) {
        float pts[20] = {0};
        int count = 0;

        float pt_vec[2][2] = {{bounds[0][0] - 1, i}, {bounds[1][0] + 1, i}};
        for (int j = 0; j < 3; j++)
        {
            float line[2][2] = {{coords[j][0],           coords[j][1]          },
                                {coords[(j + 1) % 3][0], coords[(j + 1) % 3][1]}};
            float x, y;
            if (intersection(pt_vec, line, &x, &y))
                pts[count++] = x;
        }
        if (count) {
            qsort(pts, count, sizeof(float), cmp);
            for (int j = 0; j < count; j++) {
                if (pts[j+1] && in_poly(((pts[j] + pts[j+1]) / 2), i, coords, bounds))
                    line(pts[j], i, pts[j + 1], i);  
            }
        }
    }
}

int backfaceCulling(float coords[][3], int tris[], float *visible) {
    int faces_n = 0;
    for (int i = 0; i < TRIS_SIZE; i += 3) {
        float *a = coords[tris[i]];
        float *b = coords[tris[i + 1]];
        float *c = coords[tris[i + 2]];
        float d1[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        float d2[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        float n[3] = {0};
        cross(d1, d2, n);
        if (-dot(a, n) < 0)
            visible[faces_n++] = i;
    }
    return faces_n;
}

int main() {
    float offset[3] = {0, 0, 0}, angles[3] = {0, 0, 0};
    bool too_close = false;
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
            else if (ch == 'q') {
                if (!too_close)
                    offset[2] -= speed;
            } 
            else if (ch == '6') 
                angles[2] -= speed + 10;
            else if (ch == '4') 
                angles[2] += speed + 10;
            else if (ch == '8') 
                angles[0] -= speed + 10;
            else if (ch == '2') 
                angles[0] += speed + 10;
            else if (ch == '9') 
                angles[1] -= speed + 10;
            else if (ch == '7') 
                angles[1] += speed + 10;
            else if (ch == 'x') {
                prev_coef = coef;
                if (!too_close)
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
        too_close = false;
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
        float vis_faces[FACES_SIZE] = {0};
        PerspectiveFOV(120, getwindowwidth()/getwindowheight(), 0.1, 100, projection_mat);
        float projected_coords[COORDS_SIZE][3];
        for (int i = 0; i < COORDS_SIZE; i++) {
            float projected_pt[4] = {0};
            translate(coords[i], offset);
            if (coords[i][2] < speed + 0.5)
                too_close = true;
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
        int faces_n = backfaceCulling(coords, tris, vis_faces);
        for (int i = 0; i < faces_n; i++) {
            int face_i = vis_faces[i];
            float tri[3][2] = {{projected_coords[tris[face_i]][0], projected_coords[tris[face_i]][1]},
                               {projected_coords[tris[face_i + 1]][0], projected_coords[tris[face_i + 1]][1]},
                               {projected_coords[tris[face_i + 2]][0], projected_coords[tris[face_i + 2]][1]}};
            setcolor(clrs[face_i / 3]);
            fill(tri);
        }
        setcolor(WHITE);
        // drawMesh(projected_coords, tris);
        swapbuffers();

        prev_coef = coef;
        memset(offset, 0, sizeof(offset));
        memset(angles, 0, sizeof(angles));
    }
    return 0;
}