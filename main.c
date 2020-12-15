#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

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

int add_atom(struct Atom *a, float x, float y, float z, char name[2]) {
	a->v.x = x;
	a->v.y = y;
	a->v.z = z;
	a->n_bonds = 0;
	a->lim_bonds = INITIAL_BONDS;
	a->rotate = FALSE;
	a->check = FALSE;
	strcpy(a->name, name);
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

float magnitude(struct Vector *v) {
	float n;
	n = powf(v->x, 2.) + powf(v->y, 2.) + powf(v->z, 2.);
	n = sqrtf(n);
	return n;
}

void normalise(struct Vector *v) {
	float n = magnitude(v);
	
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
	

void bond_xyz(struct Molecule *m, float bl) {
	int i, j;
	struct Vector temp;
	for (i = 0; i < m->n_atoms; i++) {
		for (j = i+1; j < m->n_atoms; j++) {
			if (i == j)
				continue;
			sub_vector(&(m->as[i].v), &(m->as[j].v), &temp);
			if (magnitude(&temp) < bl) {
				add_bond(&(m->as[i]), &(m->as[j]));
				printf("Bonding %d - %d\n", i, j);
			}
			//printf("%d - %d %f\n", i, j, magnitude(&temp));
		}
	}
	return;
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
	
	int c_line = 0;
	int c_atom = 0;
	struct Atom *ca;
	while (fgets(line, len, fp)) {
		if (c_line == 0) {
			int k = sscanf(line, "%d\n", &(m->n_atoms));
			if (k != 1) {
				fclose(fp);
				return -1;
			}
			m->as = (struct Atom *) malloc(sizeof(struct Atom) * m->n_atoms);
		} else if (c_line >= 2) {
			ca = &(m->as[c_line - 2]);
			char *token = strtok(line, " \t");
			//int k = sscanf(line, "%2s\t%f\t%f\t%f\n", ca->name, &(ca->v.x), &(ca->v.y), &(ca->v.z));
			
			//C          0.87880        0.05090       -0.00750
			int i = 0;
			float x, y, z;
			char name[2];
			while (token) {
				switch (i) {
					case 0: sprintf(name, "%2s", token); break;
					case 1: x = atof(token); break;
					case 2: y = atof(token); break;
					case 3: z = atof(token); break;
					default: fclose(fp); return -1; break;
				}
				i++;

				token = strtok(NULL, " \t");
			}
			add_atom(ca, x, y, z, name);
			printf("%s :: (%f, %f, %f) %d\n", ca->name, ca->v.x, ca->v.y, ca->v.z, c_line-2);
		}
		c_line++;
	}
	fclose(fp);
	return 1;
}

int save_xyz(struct Molecule *m, char *filename) {
	FILE * fp;
	fp = fopen(filename, "w");
	if (fp == NULL)
		return -1;
	
	fprintf(fp, "%d\n\n", m->n_atoms);
	int i;
	for (i = 0; i < m->n_atoms; i++) {
		fprintf(fp, "%2s\t%-2.6f\t%-2.6f\t%-2.6f\n", m->as[i].name, \
			m->as[i].v.x,\
			m->as[i].v.y,\
			m->as[i].v.z);
	}
	
	fclose(fp);
}

int print_dir(struct Molecule *m) {
	int i;
	for (i = 0; i < m->n_atoms; i++) {
		printf("%d\t%2s\t%-2.6f\t%-2.6f\t%-2.6f\n", i, m->as[i].name, \
			m->as[i].v.x,\
			m->as[i].v.y,\
			m->as[i].v.z);
	}
	return 0;
}

int run_script(char *filename, struct Molecule *m) {
	FILE * fp;
	char line[255];
	size_t len=255;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	
	int c_line = 0;
	int c_atom = 0;
	struct Atom *ca;
	float bl;
	int A, B;
	char command[255];
	int c;
	while (fgets(line, len, fp)) {
		reset_check(m);

		sscanf(line, "%255s %n", command, &c);
		if (strcmp(command, "open") == 0) {
			if (sscanf(line+c, "%255s", command) != 1) {
				printf("error reading script\n%s", line);
				break;
			}
			read_xyz(m, command);
		} else if (strcmp(command, "bond") == 0) {
			if (sscanf(line+c, "%f", &bl) != 1) {
				printf("Error reading script\n%s", line);
				break;
			}
			bond_xyz(m, bl);
		} else if (strcmp(command, "rotate") == 0) {
			if (sscanf(line+c, "%d %d %f", &A, &B, &bl) != 3) {
				printf("Error reading script\n%s", line);
				break;
			}
			reset_check(m);
			rotate_about(&(m->as[A]), &(m->as[B]), bl);
		} else if (strcmp(command, "output") == 0) {
			if (sscanf(line+c, "%255s", command) != 1) {
				printf("Error reading script\n%s", line);
				break;
			}
			save_xyz(m, command);
		}
	}
	fclose(fp);
	return 1;
}

int main(int argc, char *argv[]) {
	int i = 0;
	float theta = 0;

	struct Molecule mol;
	mol.n_atoms = 0;
	char command[255];
	int n, A, B;

	char script_name[255];
	if (argc >= 2) {
		strcpy(script_name, argv[1]);
		run_script(script_name, &mol);
	} else {
		do {
			reset_check(&mol);
			printf("> ");
			scanf("%255s", command);
			printf("%s\n", command);
			if (strcmp(command, "open") == 0) {
				printf(" filename > ");
				scanf("%255s", command);
				read_xyz(&mol, command);
			} else if (strcmp(command, "bond") == 0) {
				printf(" max bond length (ang) > ");
				scanf("%255s", command);
				bond_xyz(&mol, atof(command));
			} else if (strcmp(command, "graph") == 0) {
				printf(" start atom > ");
				if (scanf("%d", &n) == 1)
					print_molecule(&(mol.as[n]));
			} else if (strcmp(command, "print") == 0) {
				print_dir(&mol);
			} else if (strcmp(command, "rotate") == 0) {
				printf("This will perform a rotation about axis defined A -> B\n");
				printf("e.g. keeping A side atoms fixed, and rotating B connected atoms\n");
				printf("by an amount theta.\n");
				printf(" A > ");
				if (scanf("%s", command) != 1)
					continue;
				A = atoi(command);
				printf(" B > ");
				if (scanf("%s", command) != 1)
					continue;
				B = atoi(command);
				printf(" theta > ");
				if (scanf("%s", command) != 1)
					continue;
				theta = atof(command);
				reset_check(&mol);
				rotate_about(&(mol.as[A]), &(mol.as[B]), theta);
			} else if (strcmp(command, "output") == 0) {
				printf(" filename > ");
				if (scanf("%255s", command) != 1)
					continue;
				save_xyz(&mol, command);
			} else if (strcmp(command, "run") == 0) {
				printf(" script > ");
				if (scanf("%255s", command) != 1)
					continue;
				run_script(command, &mol);
			}
		} while (strcmp(command, "exit") != 0);
	}
	
	free_atoms(&mol);
	free(mol.as);
	return 1;
}
