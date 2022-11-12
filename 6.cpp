#define _USE_MATH_DEFINES
#include <graphics.h>
#include <math.h>
#include <float.h>

#define COORDS_SIZE 4
#define INIT_SPEED 1
#define OBJECTS_SIZE 2
#define FACES_SIZE 12
#define FACES_NUM FACES_SIZE / 3
#define CLRS_SIZE FACES_SIZE / 3
#define FACE_PARAMS 4
#define EPS 0.00001

enum face_location {INFRONT, WITHIN, BEHIND, SPLIT};
float deg2rad = M_PI / 180.0;
float ground = 1.5;

float objects[OBJECTS_SIZE][COORDS_SIZE][3] = {{{0.20, 0.5, 1.5},
                                                {0.00, 2, 2.0},
                                                {0.50, 2, 1.0},
                                                {1.00, 2, 2.0}},
    
                                               {{0.20, 0.5, 1.5},
                                                {0.00, 2, 2.0},
                                                {0.50, 2, 1.0},
                                                {1.00, 2, 2.0}}};

int graphs[FACES_SIZE] = {3, 2, 0,
                          2, 1, 0, 
                          0, 1, 3,
                          1, 2, 3};

int clrs[] = {CYAN, GREEN, RED, MAGENTA};

typedef struct tpoint {
    float x;
    float y;
    float z;
} point;

typedef struct tplane {
    float a;
    float b;
    float c;
    float d;
} plane;

typedef struct tface {
    point a;
    point b;
    point c;
    int clr;
} face;

typedef struct t_list_node {
    face *f;
    struct t_list_node *next;
} list_node;

list_node *newListNode(face tmp_f) {
    list_node *node = (list_node*)malloc(sizeof(list_node));
    face *f = (face *)malloc(sizeof(face));
    f->a = {tmp_f.a.x, tmp_f.a.y, tmp_f.a.z};
    f->b = {tmp_f.b.x, tmp_f.b.y, tmp_f.b.z};
    f->c = {tmp_f.c.x, tmp_f.c.y, tmp_f.c.z};
    f->clr = tmp_f.clr;
    node->f = f;
    node->next = NULL;
    return node;
}

typedef struct {
    list_node *head;
    list_node *tail;
} list;

list *initList() {
    list *l = (list*)malloc(sizeof(list));
	l->head = NULL;
	l->tail = NULL;
    return l;
}

int isEmpty(list *l) {
    return l->head == NULL;
}

void clearList(list *list) {
    if (list->head != NULL) {
        list_node *tmp = list->head;
        list_node *next_node;
        do {
            next_node = tmp->next;
            free(tmp->f);
            free(tmp);
        } while (tmp->next != NULL, tmp = next_node);
    }
}

void push(list *list, list_node *node) {
    if (list->head == NULL) {
        list->head = node;
        list->tail = node;
    } else {
        list->tail->next = node;
        list->tail = node;
    }
}

face pop(list *list) {
    if (list->head == NULL) {
        face res = {-1};
        return res;
    } else {
        list_node *tmp = list->head;
        face res = {0};
        res.a = {tmp->f->a.x, tmp->f->a.y, tmp->f->a.z};
        res.b = {tmp->f->b.x, tmp->f->b.y, tmp->f->b.z};
        res.c = {tmp->f->c.x, tmp->f->c.y, tmp->f->c.z};
        res.clr = tmp->f->clr;
        if (list->head == list->tail) {
            list->head = NULL;
            list->tail = NULL;
        } else 
            list->head = list->head->next;
        
        free(tmp->f);
        free(tmp);
        return res;
    }
}

typedef struct t_tree_node {
    list *faces;
    struct t_tree_node *left;
    struct t_tree_node *right;
} tree_node;

tree_node *newNode() {
    tree_node *node = (tree_node*)malloc(sizeof(tree_node));
    list *faces = initList();
    node->faces = faces;
    node->left = NULL;
    node->right = NULL;
    return node;
}

void clearTree(tree_node *root) {
    if (root == NULL) return;
    clearTree(root->left);
    clearTree(root->right);
    clearList(root->faces);
    free(root);
}

int cmp(const void * a, const void * b) {
   return ( *(float*)a - *(float*)b );
}

bool ccw(float *a, float *b, float *c) {
    return (c[1]-a[1])*(b[0]-a[0]) > (b[1]-a[1])*(c[0]-a[0]);
}

bool intersect(float a[2], float b[2], float c[2], float d[2]) {
    return (ccw(a, c, d) != ccw(b, c, d)) && (ccw(a, b, c) != ccw(a, b, d));
}

float dot(point a, point b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

point cross(point a, point b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

point sum(point a, point b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

point subtract(point a, point b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

point mult(point a, float n) {
    return {a.x*n, a.y*n, a.z*n};
}

point divide(point a, float n) {
    return {a.x / n, a.y / n, a.z / n};
}

bool intersection(float line1[2][2], float line2[2][2], float *i_x, float *i_y) {
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

plane face2plane(face f) {
    point d1 = subtract(f.b, f.a);
    point d2 = subtract(f.c, f.a);
    point n = cross(d1, d2);
    return {n.x, n.y, n.z, -dot(f.a, n)};
}
plane face2plane(face *f) {
    point d1 = subtract(f->b, f->a);
    point d2 = subtract(f->c, f->a);
    point n = cross(d1, d2);
    return {n.x, n.y, n.z, -dot(f->a, n)};
}
float dist(point pt, plane p) {
    return p.a*pt.x + p.b*pt.y + p.c*pt.z + p.d;
}

enum face_location located(plane p, face f) {
    int has_pt_infront = 0, has_pt_behind = 0;
    point poly[3] = {f.a, f.b, f.c};
    for (int i = 0; i < 3; i++) {
        float sign = dist(poly[i], p);
        if (sign > EPS) has_pt_infront = 1; 
        if (sign < -EPS) has_pt_behind = 1;
    }
    if(has_pt_infront && has_pt_behind) return SPLIT;
    if(has_pt_infront) return INFRONT;
    if(has_pt_behind) return BEHIND;
    return WITHIN;
}
//unused
void planesIntersection(plane p1, plane p2, point *r_pt, point *r_dir) {
    point n1 = {p1.a, p1.b, p1.c}, n2 = {p2.a, p2.b, p2.c};
    point n3 = cross(n1, n2);
    float det = n3.x*n3.x + n3.y*n3.y + n3.z*n3.z;
    point tmp_pt = divide(sum(mult(cross(n3, n2), p1.d),
                              mult(cross(n1, n3), p2.d)), det);
    r_pt->x = tmp_pt.x; r_dir->x = n3.x;
    r_pt->y = tmp_pt.y; r_dir->y = n3.y;
    r_pt->z = tmp_pt.z; r_dir->z = n3.z;
}

point isect(point start, point end, plane p) {
    point ray = subtract(end, start);
    point n = {p.a, p.b, p.c}; 
    if (dot(n, ray) == 0)
        return {FLT_MAX, FLT_MAX, FLT_MAX};
    float t = -(p.d +dot(start, n)) / dot(n, ray);
    return {start.x + t*ray.x, start.y + t*ray.y, start.z + t*ray.z};
}

void split(face cur_f, face f, list *unsorted_faces) {
    plane p = face2plane(cur_f);

    point poly[3] = {f.a, f.b, f.c};
    point i[3] = {0}, b[3] = {0};
    int i_num = 0, b_num = 0;
    for (int k = 0; k < 3; k++) {
        float d = dist(poly[k], p);
        if (d > 0) i[i_num++] = poly[k];
        else b[b_num++] = poly[k];
    }
    
    if (i_num == 1) {
        face new_faces[] = {{i[0], isect(i[0], b[0], p), isect(i[0], b[1], p), f.clr},
                            {b[0], b[1], isect(b[0], i[0], p), f.clr},
                            {b[1], isect(b[0], i[0], p), isect(b[1], i[0], p), f.clr}};
        for (int i = 0; i < 3; i++)
            push(unsorted_faces, newListNode(new_faces[i]));
    } else if (i_num == 2) {
        face new_faces[] = {{b[0], isect(b[0], i[0], p), isect(b[0], i[1], p), f.clr},
                            {i[0], i[1], isect(i[0], b[0], p), f.clr},
                            {i[1], isect(i[0], b[0], p), isect(i[1], b[0], p), f.clr}};
        for (int i = 0; i < 3; i++)
            push(unsorted_faces, newListNode(new_faces[i]));
    } 
}

void bspSort(tree_node *root, list *unsorted_faces) {
    list *left = initList(), *current = initList(), *right = initList();
    face cur_f = pop(unsorted_faces);
    push(current, newListNode(cur_f));
    while (!isEmpty(unsorted_faces)) {
        face f = pop(unsorted_faces);
        enum face_location loc = located(face2plane(cur_f), f);
        if (loc == INFRONT)
            push(left, newListNode(f));
        else if (loc == WITHIN)
            push(current, newListNode(f));
        else if (loc == BEHIND)
            push(right, newListNode(f));
        else if (loc == SPLIT) {
            split(cur_f, f, unsorted_faces);
        }
    }
    root->faces = current;
    if (!isEmpty(left)) {
        root->left = newNode();
        bspSort(root->left, left);
    } else 
        root->left = NULL;
    if (!isEmpty(right)) {
        root->right = newNode();
        bspSort(root->right, right);
    } else
        root->right = NULL;
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

bool in_poly(float x, float y, float coords[3][3], float bounds[2][2]) {
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

void fill(float coords[3][3], int clr) {
    setcolor(clr);
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
    setcolor(WHITE);
}

void project(point *poly, float projected_coords[3][3]) {
    float projection_mat[16] = {0};
    PerspectiveFOV(90, getwindowwidth()/getwindowheight(), 0.1, 50, projection_mat);
    for (int i = 0; i < 3; i++) {
        float pt[] = {poly[i].x, poly[i].y, poly[i].z, 0};
        float projected_pt[4] = {0};
        mult(pt, projection_mat, projected_pt);
        for (int j = 0; j < 3; j++)
            projected_coords[i][j] = 0.5 * (projected_pt[j] / projected_pt[3] + 1);
        projected_coords[i][0] *= getwindowwidth();
        projected_coords[i][1] *= getwindowheight();
    }
}

void render(face f) {
    point poly[3] = {f.a, f.b, f.c};
    float projected_coords[3][3];
    project(poly, projected_coords);
    fill(projected_coords, f.clr);

    // for (int i = 0; i < 3; i++) {
    //     float *pt1 = projected_coords[i];
    //     float *pt2 = projected_coords[(i + 1) % 3];
    //     line(pt1[0], pt1[1], pt2[0], pt2[1]);
    // }
}

void paint(tree_node *root, point view) {
    if(root == NULL)
        return;
    plane p = face2plane(root->faces->head->f);
    float view_sign = dist(view, p);
    if(view_sign > EPS){
        paint(root->right, view);
        while(!isEmpty(root->faces))
            render(pop(root->faces));
        paint(root->left, view);
    } else {
        paint(root->left, view);
        while(!isEmpty(root->faces))
            render(pop(root->faces));
        paint(root->right, view);
    }
    
}

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

void drawMesh(float coords[][3], int graphs[]) {
    for (int i = 0; i < FACES_SIZE; i += 3) {
        for (int j = 0; j < 3; j++) {
            float *pt1 = coords[graphs[i + j]];
            float *pt2 = coords[graphs[i + (j + 1) % 3]];
            line(pt1[0], pt1[1], pt2[0], pt2[1]);
        }
    }
}

point convert(float pt[3]) {
    point res = {pt[0], pt[1], pt[2]};
    return res;
}

void bakeShadow(float coords[][3], int graphs[], point light_pt, list *visible) {
    for (int i = 0; i < FACES_SIZE; i += 3) {
        point a = convert(coords[graphs[i]]);
        point b = convert(coords[graphs[i + 1]]);
        point c = convert(coords[graphs[i + 2]]);
        point poly[3] = {a, b, c};
        point shadow[3];
        float projected_coords[3][3];
        for (int i = 0; i < 3; i++) {
            float t = (ground - poly[i].y) / light_pt.y;
            shadow[i] = {poly[i].x + t*(poly[i].x - light_pt.x), 
                         ground,
                         poly[i].z + t*(poly[i].x - light_pt.z)};
        }
        push(visible, newListNode({shadow[0], shadow[1], shadow[2], LIGHTGRAY}));
    }
}

void backfaceCulling(float coords[][3], int graphs[], list *visible) {
    for (int i = 0; i < FACES_SIZE; i += 3) {
        point a = convert(coords[graphs[i]]);
        point b = convert(coords[graphs[i + 1]]);
        point c = convert(coords[graphs[i + 2]]);
        point d1 = subtract(b, a);
        point d2 = subtract(c, a);
        point n = cross(d1, d2);
        if (dot(a, n) > 0) {
            face f = {a, b, c, clrs[i / 3]};
            push(visible, newListNode(f));
        }
    }
}

int main() {
    point view = {0, 0, 0};
    point light_pt = {20, 10, 0};
    float offsets[2][3] = {{-1, -0.5, 0}, {-1, -1, 0}}; 
    float angles[2][3]  = {{0, -80, -40}, {0, 0, 0}};
    bool too_close = false;
    int angle = 0, n = 0;
    float speed = INIT_SPEED;
    float coefs[2][2] = {{1, 1}, {1, 1}};
    int gd = DETECT, gm;
    initgraph(&gd, &gm, (char *)"");
    while (1) {
        if (kbhit()) {
            char ch = getch();
            float *offset = offsets[n];
            float *angle = angles[n];
            float *coef = coefs[n];
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
            else if (ch == 'k')
                light_pt.y -= speed;
            else if (ch == 'l')
                light_pt.y += speed;
            else if (ch == 'h')
                light_pt.x -= speed;
            else if (ch == 'j')
                light_pt.x += speed;
            else if (ch == 'u')
                light_pt.z += speed;
            else if (ch == 'i') {
                light_pt.z -= speed;
            }
            else if (ch == '6') 
                angle[2] -= speed + 10;
            else if (ch == '4') 
                angle[2] += speed + 10;
            else if (ch == '8') 
                angle[0] -= speed + 10;
            else if (ch == '2') 
                angle[0] += speed + 10;
            else if (ch == '9') 
                angle[1] -= speed + 10;
            else if (ch == '7') 
                angle[1] += speed + 10;
            else if (ch == '0')
                n = 0;
            else if (ch == '1')
                n = 1;
            else if (ch == 'x') {
                coef[0] = coef[1];
                if (!too_close)
                    coef[1] += 0.01 * speed;
            }
            else if (ch == 'z') {
                coef[0] = coef[1];
                if (coef[1] > 0.01 * speed)
                    coef[1] -= 0.01 * speed;
            }
            else if (ch == '+')
                speed += 1;
            else if (ch == '-')
                if (speed > 1)
                    speed -= 1;
        }
        too_close = false;
        cleardevice();
        
        list *unsorted_faces = initList();
        float (*objects_p)[COORDS_SIZE][3] = objects;
        for (int k = 0; k < OBJECTS_SIZE; k++) {
            float (*coords)[3] = objects_p[k];
            float *offset = offsets[k];
            float *angle = angles[k];
            float origin[3] = {0}, neg_origin[3] = {0};
            for (int i = 0; i < COORDS_SIZE; i++)
                for (int j = 0; j < 3; j++)
                    origin[j] += coords[i][j];
            for (int i = 0; i < 3; i++) {
                origin[i] /= COORDS_SIZE;
                neg_origin[i] = -origin[i];
            }

            for (int i = 0; i < COORDS_SIZE; i++) {
                float projected_pt[4] = {0};
                translate(coords[i], offset);
                if (coords[i][2] < speed + 0.5 && k == n)
                    too_close = true;
                translate(coords[i], neg_origin);
                rotate(coords[i], angle);
                scale(coords[i], coefs[k][0], coefs[k][1]);
                translate(coords[i], origin);
            }
            bakeShadow(coords, graphs, light_pt, unsorted_faces);
            backfaceCulling(coords, graphs, unsorted_faces);
        }
        tree_node *root = newNode();
        bspSort(root, unsorted_faces);
        paint(root, view);

        swapbuffers();
        for (int i = 0; i < OBJECTS_SIZE; i++) {
            coefs[i][0] = coefs[i][1];
            memset(offsets[i], 0, sizeof(offsets[i]));
            memset(angles[i], 0, sizeof(angles[i]));
        }
    }
    return 0;
}