/*
Author: Seoyun Lee						Class: CPSC 445
Homework Assignment: 3

uses a quad tree to find the k nearest neighbors of a point.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

typedef struct point {
	double x;
	double y;
	int index;
} point;

typedef struct node {
	struct node *previous;
	struct node *t_left;
	struct node *t_right;
	struct node *b_left;
	struct node *b_right;
	int num_points;
	double side_length;
	point *t_left_point;
	point *t_right_point;
	point *b_left_point;
	point *b_right_point;
	point **points;
} node;

/* declare functions */
point *create_point(double x, double y, int index);
node *create_node(node *parent, int number, point *b_left, point **points_contained, double length);
double distance(point *a, point *b);
double find_max(double a, double b);
void expand_tree(node *parent, int k);
int which_box(node *parent, point *current_point);
int intersects(point *a, double radius, node *current);
int is_in_box(double a_x, double a_y, node *current);
void search_for_points(node *box, double *distances, int k, point *a, double radius, int *iz);
void search_tree(node *current, int k, int *iz);
void find_nearest_points(node *leaf, point *a, int k, int *iz);
void add_point(node *box, double *distances, int k, point *a, int *iz);
void seek(double *a, int n, int k, int *iz);
void seek_naive(double *a, int n, int k, int *iz);

int main() {
    return 0;
}

/* creates a point */
point *create_point(double x, double y, int index) {
	point *new_pt = malloc(sizeof(point));
	new_pt->x = x;
	new_pt->y = y;
	new_pt->index = index;
	return new_pt;
}

/* creates a node */
node *create_node(node *parent, int number, point *b_left, point **points_contained, double length) {
	node *new_node = malloc(sizeof(node));
	new_node->previous = parent;
	new_node->t_left = NULL;
	new_node->t_right = NULL;
	new_node->b_left = NULL;
	new_node->b_right = NULL;
	new_node->num_points = number;
	new_node->side_length = length;
	new_node->t_left_point = create_point((*b_left).x, (*b_left).y + length, -1);
	new_node->t_right_point = create_point((*b_left).x + length, (*b_left).y + length, -1);
	new_node->b_left_point = create_point((*b_left).x, (*b_left).y, -1);
	new_node->b_right_point = create_point((*b_left).x + length, (*b_left).y, -1);
	new_node->points = points_contained;
	return new_node;
}

/*
	distance between two points a and b
*/
double distance(point *a, point *b) {
	double dx, dy;

	dx = b->x - a->x;
	dy = b->y - a->y;

	return sqrt(dx * dx + dy * dy);
}

double find_max(double a, double b) {
	if(a > b) {
		return a;
	}
	else {
		return b;
	}
}

/* expands out the tree recursivey */
void expand_tree(node *parent, int k) {
	/* stop if there are less than k points */
	if(parent->num_points <= k) {
		return;
	}

	point **b_left_points = malloc(parent->num_points * sizeof(point *));
	point **b_right_points = malloc(parent->num_points * sizeof(point *));
	point **t_left_points = malloc(parent->num_points * sizeof(point *));
	point **t_right_points = malloc(parent->num_points * sizeof(point *));

	int size_b_left = 0, size_b_right = 0, size_t_left = 0, size_t_right = 0;

	int select_box, i;
	for(i = 0; i < parent->num_points; i++) {
		select_box = which_box(parent, parent->points[i]);
		if(select_box == 0) {
			b_left_points[size_b_left] = parent->points[i];
			size_b_left++;
		}
		else if(select_box == 1) {
			b_right_points[size_b_right] = parent->points[i];
			size_b_right++;
		}
		else if(select_box == 2) {
			t_right_points[size_t_right] = parent->points[i];
			size_t_right++;
		}
		else {
			t_left_points[size_t_left] = parent->points[i];
			size_t_left++;
		}
	}

	double new_length;
	new_length = (parent->side_length)/2.0;

	parent->b_left = create_node(parent, size_b_left, parent->b_left_point, b_left_points, new_length);
	parent->b_right = create_node(parent, size_b_right, create_point(parent->b_left_point->x + new_length, parent->b_left_point->y, -1), b_right_points, new_length);
	parent->t_right = create_node(parent, size_t_right, create_point(parent->b_left_point->x + new_length, parent->b_left_point->y + new_length, -1), t_right_points, new_length);
	parent->t_left = create_node(parent, size_t_left, create_point(parent->b_left_point->x, parent->b_left_point->y + new_length, -1), t_left_points, new_length);

	expand_tree(parent->b_left, k);
	expand_tree(parent->b_right, k);
	expand_tree(parent->t_right, k);
	expand_tree(parent->t_left, k);
}

/* 0 is bottom left, 1 is bottom right, 2 is top right, 3 is top left */
int which_box(node *parent, point *current_point) {
	double new_length;
	new_length = (parent->side_length) / 2.0;
	if(current_point->x < (parent->b_left_point->x + new_length) 
		&& current_point->y < (parent->b_left_point->y + new_length)) {
	   	return 0;
	}
	else if(current_point->x >= (parent->b_left_point->x + new_length)
		&& current_point->y < (parent->b_left_point->y + new_length)) {
		return 1;
	}
	else if(current_point->x >= (parent->b_left_point->x + new_length)
		&& current_point->y >= (parent->b_left_point->y + new_length)) {
		return 2;
	}
	else {
		return 3;
	}
}

/* seaches the tree for leaves and then finds the nearest points for all points in each leaf */
void search_tree(node *current, int k, int *iz) {
	/* we are at a leaf */
	if(current->t_left == NULL) {
		int i;
		for(i = 0; i < current->num_points; i++) {
			find_nearest_points(current, current->points[i], k, iz);
		}
		return;
	}
	/* not at a leaf, keep going until we find one */
	search_tree(current->t_left, k, iz);
	search_tree(current->t_right, k, iz);
	search_tree(current->b_left, k, iz);
	search_tree(current->b_right, k, iz);
}

/* finds the k nearest points to point a */
void find_nearest_points(node *leaf, point *a, int k, int *iz) {
	double radius;
	radius = distance(a, leaf->previous->t_left_point);
	radius = find_max(radius, distance(a, leaf->previous->t_right_point));
	radius = find_max(radius, distance(a, leaf->previous->b_left_point));
	radius = find_max(radius, distance(a, leaf->previous->b_right_point));

	node *largest_box;
	largest_box = (node*) malloc(sizeof(node));
	largest_box = leaf->previous;
	while(intersects(a, radius, largest_box) == 1 && largest_box->previous != NULL) {
		largest_box = largest_box->previous;
	}

	/* given the largest box, search its children recursively */
	double *distances;
	distances = malloc(k * sizeof(double));
	int i;
	for(i = 0; i < k; i++) {
		distances[i] = 1.0e21;
	}
	search_for_points(largest_box, distances, k, a, radius, iz);
	free(distances);
}

/* adds point to iz */
void add_point(node *box, double *distances, int k, point *a, int *iz) {
	int i, j;
	for(i = 0; i < box->num_points; i++) {
		double dist;
		dist = distance(a, box->points[i]);
		if(dist == 0) {
			continue;
		}
		j = 0;
		if(dist < distances[j]) {
			while(dist < distances[j] && j < k-1 ) {
				if(dist < distances[j + 1]) {
					distances[j] = distances[j + 1];
					iz[(a->index)*k + j] = iz[(a->index)*k + j + 1];
					j++;
				}
				else {
					distances[j] = dist;
					break;
				}
			}
			distances[j] = dist;
			iz[(a->index)*k + j] = box->points[i]->index;
		}
	}
}

/* searches for viable points */
void search_for_points(node *box, double *distances, int k, point *a, double radius, int *iz) {
	/* check if we're at a leaf */
	if(box->t_left == NULL) {
		add_point(box, distances, k, a, iz);
		return;
	}
	if(intersects(a, radius, box->t_left)) {
		search_for_points(box->t_left, distances, k, a, radius, iz);
	}
	if(intersects(a, radius, box->t_right)) {
		search_for_points(box->t_right, distances, k, a, radius, iz);
	}
	if(intersects(a, radius, box->b_left)) {
		search_for_points(box->b_left, distances, k, a, radius, iz);
	}
	if(intersects(a, radius, box->b_right)) {
		search_for_points(box->b_right, distances, k, a, radius, iz);
	}
	else {
		return;
	}
}

/* does the circle intersect with a box? */
int intersects(point *a, double radius, node *current) {
	if(distance(a, current->t_left_point) <= radius ||
		distance(a, current->t_right_point) <= radius ||
		distance(a, current->b_left_point) <= radius ||
		distance(a, current->b_right_point) <= radius) {
		return 1;
	}
	if(is_in_box(a->x + radius, a->y, current) == 1 || is_in_box(a->x - radius, a->y, current) == 1 ||
		is_in_box(a->x, a->y + radius, current) == 1 || is_in_box(a->x, a->y - radius, current) == 1) {
		return 1;
	}
	if(is_in_box(a->x, a->y, current) == 1) {
		return 1;
	}
	else {
		return 0;
	}
}

/* decides if a point is in a box */
int is_in_box(double a_x, double a_y, node *current) {
	if(a_x <= current->t_right_point->x && a_x >= current->t_left_point->x
		&& a_y <= current->t_left_point->y && a_y >= current->b_left_point->y) {
		return 1;
	}
	else {
		return 0;
	}
}

void seek(double *a, int n, int k, int *iz) {
	/* store a as array of points */
	point **points;
	points = malloc(n * sizeof(point*));

	int i;
	for(i = 0; i < n; i++) {
		points[i] = create_point(a[2 * i], a[2 * i + 1], i);
	}

	/* find min and max points */
	double min_x = points[0]->x, max_x = points[0]->x;
	double min_y = points[0]->y, max_y = points[0]->y;

	for(i = 1; i < n; i++) {
		if(points[i]->x < min_x) {
			min_x = points[i]->x;
		}
		if(points[i]->x > max_x) {
			max_x = points[i]->x;
		}
		if(points[i]->y < min_y) {
			min_y = points[i]->y;
		}
		if(points[i]->y > max_y) {
			max_y = points[i]->y;
		}
	}

	point *b_left;
	b_left = (point*) malloc(sizeof(point));
	b_left = create_point(min_x, min_y, -1);

	double side_length;
	side_length = find_max(max_x - min_x, max_y - min_y);

	node *root;
	root = (node*) malloc(sizeof(node));
	root = create_node(NULL, n, b_left, points, side_length);

	expand_tree(root, k);
	search_tree(root, k, iz);

	free(points);
	free(b_left);
	free(root);
}


void seek_naive(double *a, int n, int k, int *iz) {
	point **points;
	points = malloc(n * sizeof(point*));

	int i, j;
	for(i = 0; i < n; i++) {
		points[i] = create_point(a[2 * i], a[2 * i + 1], i);
	}

	for(i = 0; i < n; i++) {
		double *distances;
		distances = malloc(k * sizeof(double));
		int l;
		for(l = 0; l < k; l++) {
			distances[l] = 1.0e21;
		}

		for(j = 0; j < n; j++) {
			if(i == j) {
				continue;
			}
			double dist;
			dist = distance(points[i], points[j]);
			l = 0;
			if(dist < distances[l]) {
				while(dist < distances[l] && l < k-1) {
					if(dist < distances[l + 1]) {
						distances[l] = distances[l + 1];
						iz[(points[i]->index)*k + l] = iz[(points[i]->index)*k + l + 1];
						l++;
					}
					else {
						distances[l] = dist;
						break;
					}
				}
				distances[l] = dist;
				iz[(points[i]->index)*k + l] = points[j]->index;
			}

		}
		free(distances);
	}
}




