#define _USE_MATH_DEFINES
#include <stdio.h>
#include <graphics.h>
#include <math.h>

#define COORDS_SIZE 4
#define INIT_SPEED 1
#define OBJECTS_SIZE 2
#define TRIS_SIZE 12
#define FACES_SIZE TRIS_SIZE / 3

float objects[OBJECTS_SIZE][COORDS_SIZE][3] = {{{0.20, 0.5, 1.5},
                                                {-0.00, 2, 2.0},
                                                {-0.50, 2, 1.0},
                                                {-1.00, 2, 2.0}},
    
                                               {{0.20, 0.5, 1.5},
                                                {0.00, 2, 2.0},
                                                {0.50, 2, 1.0},
                                                {1.00, 2, 2.0}}};

struct node {
    int *face;
    struct node *left;
    struct node *right;
};

struct node *newNode(int a, int b, int c) {
    struct node *node = (struct node*)malloc(sizeof(struct node));
    int *face = (int *)malloc(3*sizeof(int));
    face[0] = a; face[1] = b; face[2] = c;
    node->face = face;
    node->left = NULL;
    node->right = NULL;
    return node;
}

int main() {
    struct node *root = newNode(1, 2, 3);
    root->left = newNode(2, 2, 2);
}