#define _USE_MATH_DEFINES
#include <stdio.h>
#include <graphics.h>
#include <math.h>

#define COORDS_SIZE 4
#define INIT_SPEED 1
#define OBJECTS_SIZE 2
#define TRIS_SIZE 12
#define FACES_SIZE TRIS_SIZE / 3

typedef struct tpoint {
    float x;
    float y;
    float z;
} point;

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

point objects[OBJECTS_SIZE][COORDS_SIZE] = {{{0.20, 0.5, 1.5},
                                                {0.00, 2, 2.0},
                                                {0.50, 2, 1.0},
                                                {1.00, 2, 2.0}},
    
                                               {{0.20, 0.5, 1.5},
                                                {0.00, 2, 2.0},
                                                {0.50, 2, 1.0},
                                                {1.00, 2, 2.0}}};

int main() {
    list *l = initList();
    {
        face f = {{1, 2, 3,}, {1, 2, 3}, {1, 2, 3}, 0};
        push(l, newListNode(f));
    }
    face f = pop(l);
}