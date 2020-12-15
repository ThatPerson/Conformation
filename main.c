#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0
#define INITIAL_BONDS 5

struct Vector {
	float x, y, z;
};

struct Atom {
	struct Vector v;
	int rotate;
	int check;
	int n_bonds;
	int lim_bonds;
	struct Atom ** bonds;
	char name[2];
};



struct Molecule {
	struct Atom *as;
	int n_atoms;
};

int add_atom(struct Atom *a, float x, float y, float z) {
	a->v.x = x;
	a->v.y = y;
	a->v.z = z;
	a->n_bonds = 0;
	a->lim_bonds = INITIAL_BONDS;
	a->rotate = FALSE;
	a->check = FALSE;
	a->bonds = (struct Atom **) malloc(sizeof(struct Atom *) * INITIAL_BONDS);
	return 1;
}

void reset_check(struct Molecule *m) {
	int i;
	for (i = 0; i < m->n_atoms; i++) {
		m->as[i].check = FALSE;
	}
	return;
}

int add_bond(struct Atom *a, struct Atom *b) {
	a->bonds[a->n_bonds] = b;
	b->bonds[b->n_bonds] = a;
	a->n_bonds++;
	b->n_bonds++;
	if (a->n_bonds >= a->lim_bonds) {
		a->lim_bonds *= 2;
		a->bonds = (struct Atom **) realloc(a->bonds, sizeof(struct Atom *) * a->lim_bonds);
	}
	if (b->n_bonds >= b->lim_bonds) {
		b->lim_bonds *= 2;
		b->bonds = (struct Atom **) realloc(b->bonds, sizeof(struct Atom *) * b->lim_bonds);
	}
	return 1;
}

void free_atoms(struct Molecule *m) {
	int i;
	for (i = 0; i < m->n_atoms; i++) {
		free(m->as[i].bonds);
	}
	return;
}

void print_moleculef(struct Atom * a, int n) {
	if (a->check == TRUE)
		return;
	a->check = TRUE;
	int i,j;
	for (j = 0; j < n+1; j++) 
		printf("> ");
	printf("Atom [%f, %f, %f] %d\n", a->v.x, a->v.y, a->v.z, a->n_bonds);
	
	for (i = 0; i < a->n_bonds; i++) {
		print_moleculef(a->bonds[i], n+1);
	}
	return;
}

void print_molecule(struct Atom * a) {
	if (a->check == TRUE)
		return;
	a->check = TRUE;
	int i;
	printf("%f, %f, %f\n", a->v.x, a->v.y, a->v.z);
	
	for (i = 0; i < a->n_bonds; i++) {
		print_molecule(a->bonds[i]);
	}
	return;
}

int rodrigues_rotation(struct Vector *a, struct Vector *v, float theta) {
	struct Vector or;
	or.x = a->x;
	or.y = a->y;
	or.z = a->z;
	
	float R[3][3];
	float ct = cos(theta), st = sin(theta);
	R[0][0] = cos(theta) + powf(v->x, 2.) * (1-ct);
	R[0][1] = v->x * v->y * (1-ct) - v->z * st;
	R[0][2] = v->y * st + v->x * v->z * (1 - ct);
	R[1][0] = v->z * st + v->x * v->y * (1 - ct);
	R[1][1] = ct + powf(v->y, 2.) * (1 - ct);
	R[1][2] = -v->x * st + v->y * v->z * (1 - ct);
	R[2][0] = -v->y * st + v->x * v->z * (1 - ct);
	R[2][1] = v->x * st + v->y * v->z * (1 - ct);
	R[2][2] = ct + powf(v->z, 2.) * (1 - ct);
	
	a->x = R[0][0] * or.x + R[0][1] * or.y + R[0][2] * or.z;
	a->y = R[1][0] * or.x + R[1][1] * or.y + R[1][2] * or.z;
	a->z = R[2][0] * or.x + R[2][1] * or.y + R[2][2] * or.z;
	
	return 1;
}

void sub_vector(struct Vector *a, struct Vector *b, struct Vector *res) {
	// res = a - b
	res->x = a->x - b->x;
	res->y = a->y - b->y;
	res->z = a->z - b->z;
	return;
}

void add_vector(struct Vector *a, struct Vector *b, struct Vector *res) {
	// res = a + b
	res->x = a->x + b->x;
	res->y = a->y + b->y;
	res->z = a->z + b->z;
	return;
}

void rotate(struct Atom *a, struct Vector *v, struct Vector * offset, float theta) {
	if (a->check == TRUE)
		return;
	a->check = TRUE;
	sub_vector(&(a->v), offset, &(a->v));
	rodrigues_rotation(&(a->v), v, theta);
	add_vector(&(a->v), offset, &(a->v));
	int i;
	for (i = 0; i < a->n_bonds; i++) {
		rotate(a->bonds[i], v, offset, theta);
	}
	return;
}

void normalise(struct Vector *v) {
	float n;
	n = powf(v->x, 2.) + powf(v->y, 2.) + powf(v->z, 2.);
	n = sqrtf(n);
	v->x /= n;
	v->y /= n;
	v->z /= n;
	return;
}

int rotate_about(struct Atom *a, struct Atom *b, float theta) {
	// rotate about axis a->b, varying b neighbours.
	// need to perform reset check before running this code or BAD THINGS
	// may happen.
	
	a->check = TRUE;
	a->rotate = FALSE;
	
	struct Vector offset;
	offset.x = a->v.x;
	offset.y = a->v.y;
	offset.z = a->v.z;
	
	struct Vector rotate_about;
	sub_vector(&(b->v), &(a->v), &rotate_about);
	normalise(&rotate_about);
	rotate(b, &rotate_about, &offset, theta);
	return 1;
}

int read_xyz(struct Molecule *m, char *filename) {
	FILE * fp;
	char line[255];
	size_t len=255;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}

	while (fgets(line, len, fp)) {
		printf("%s\n", line);
	}
	fclose(fp);
	return 1;
}

int main(int argc, char *argv[]) {
	int i = 0;
	float theta = 0;
	if (argc >= 2)
		i = atoi(argv[1]);
	if (argc >= 3)
		theta = atof(argv[2]);
		
	struct Molecule mol;
	
	//read_xyz(&mol, "mp.xyz");
	
	//return -1;
	mol.n_atoms = 4; // as illustrative example
	mol.as = (struct Atom *) malloc(sizeof(struct Atom) * mol.n_atoms);
	add_atom(&(mol.as[0]), 0, 0, 0);
	add_atom(&(mol.as[1]), 1, 1, 0);
	add_atom(&(mol.as[2]), 1, 2, 0);
	add_atom(&(mol.as[3]), 2, 2, 0);
	add_bond(&(mol.as[0]), &(mol.as[1]));
	add_bond(&(mol.as[1]), &(mol.as[2]));
	add_bond(&(mol.as[2]), &(mol.as[3]));
	
	printf("Pre Rotation\n");
	print_molecule(&(mol.as[i]));
	reset_check(&mol);
	rotate_about(&(mol.as[1]), &(mol.as[2]), theta);
	reset_check(&mol);
	printf("Post Rotation\n");
	print_molecule(&(mol.as[i]));
	free_atoms(&mol);
	free(mol.as);
	return 1;
}

