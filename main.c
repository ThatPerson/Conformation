#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define TRUE 1
#define FALSE 0
#define INITIAL_BONDS 5 // Number of bonds per atom possible.
						// If you've got some weird compound may need to increase

struct Vector {
	float x, y, z;
};

/* v is position, rotate and check are TRUE or FALSE and should be
   reset post processing. name is atom name, bonds is pointer array to other atoms */
struct Atom {
	struct Vector v;
	int rotate;
	int check;
	int n_bonds;
	int lim_bonds;
	struct Atom ** bonds;
	char name[2];
};

/* Molecule container - n_atoms in as struct. */
struct Molecule {
	struct Atom *as;
	int n_atoms;
};

/* add_atom()
 *  Initialises atom *a.
 * arg *a: Pointer to atom struct to put data in.
 * arg x, y, z: positional coordinated.
 * arg name: atom name
 * returns: 1 in all cases.
 
 * Error modes;
 *  - If atom at *a not allocated memory, will segfault.
 */
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

/* reset_check()
 *  Resets all .check values to be FALSE. Should be run before _any_
 *  graph or rotate commands
 
 * Error modes;
 *  - If not all atoms are allocated, will segfault.
 */
void reset_check(struct Molecule *m) {
	int i;
	for (i = 0; i < m->n_atoms; i++) {
		m->as[i].check = FALSE;
	}
	return;
}

/* add_bond()
 *  Creates bond between atom *a and *b. 
 * arg *a: Pointer to atom.
 * arg *b: Pointer to atom.
 * returns: 1 in all cases
 
 * Error modes;
 *  - Will seg fault if *a or *b do not exist, or if there is insufficient memory.
 */
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

/* free_atoms()
 *  Frees bond arrays for each atom within molecule *m
 *
 * Error modes;
 *  - Will double free if any atom has not been initialized.
 */
void free_atoms(struct Molecule *m) {
	int i;
	for (i = 0; i < m->n_atoms; i++) {
		free(m->as[i].bonds);
	}
	return;
}

/* print_moleculef()
 *  Prints molecular graph starting at atom n.
 *  Iterates in a functional manner, setting check to TRUE once an
 *  atom has been visited. Visits each atom once.
 * arg *a: Pointer to initial atom.
 * arg  n: Atom to begin graph at
 */
void print_moleculef(struct Atom * a, int n) {
	// if we have been visited, terminate execution.
	if (a->check == TRUE)
		return;
	a->check = TRUE;
	int i,j;
	for (j = 0; j < n+1; j++) 
		printf("> ");
	printf("Atom %s [%f, %f, %f] %d\n", a->name, a->v.x, a->v.y, a->v.z, a->n_bonds);
	
	// Recursively loop over, calling self for each neighbour.
	for (i = 0; i < a->n_bonds; i++) {
		print_moleculef(a->bonds[i], n+1);
	}
	return;
}

/* print_molecule()
 *  Prints molecular system, starting with one atom and looping to others
 */
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

/* rodrigues_rotation()
 *  Applies rotation of theta about vector *v to vector *a (in situ).
 *  Vector *v should be from offset (eg if defined as bond, you need to
 *  subtract that offset from vector *a prior to applying rotation).
 * arg *a: Vector to be rotated.
 * arg *v: Vector to rotate about.
 * arg theta: Angle (in radians)
 * returns: 1 in all cases
 */
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

/* sub_vector()
 *  Subtracts vector *b from *a, and outputs into *res
 * arg *a: Vector
 * arg *b: Vector
 * arg *res: Vector, res = a - b
 */
void sub_vector(struct Vector *a, struct Vector *b, struct Vector *res) {
	// res = a - b
	res->x = a->x - b->x;
	res->y = a->y - b->y;
	res->z = a->z - b->z;
	return;
}

/* add_vector()
 *  Adds vector *b to *a, and outputs into *res
 * arg *a: Vector
 * arg *b: Vector
 * arg *res: Vector, res = a + b
 */
void add_vector(struct Vector *a, struct Vector *b, struct Vector *res) {
	// res = a + b
	res->x = a->x + b->x;
	res->y = a->y + b->y;
	res->z = a->z + b->z;
	return;
}

/* rotate()
 *  Provides interface to rodrigues_rotation() with offset. 
 *  Rotates atom *a about vector *v an amount theta radians, taking into
 *  account an offset *offset. Will then loop over all atoms bonded to atom *a
 *  and rotate these too, unless their check flag has been enabled (e.g. they have
 *  been rotationally unallowed or have already been rotated.)
 * arg *a: Atom to be rotated.
 * arg *v: Vector about which to rotate
 * arg *offset: Vector to offset by.
 * arg theta: Angle to rotate by in radians
 */
void rotate(struct Atom *a, struct Vector *v, struct Vector * offset, float theta) {
	if (a->check == TRUE)
		return;
	a->check = TRUE;
	// we subtract the offset away, 
	sub_vector(&(a->v), offset, &(a->v));
	// then perform the rotation,
	rodrigues_rotation(&(a->v), v, theta);
	// then add the offset back on
	add_vector(&(a->v), offset, &(a->v));
	int i;
	for (i = 0; i < a->n_bonds; i++) {
		rotate(a->bonds[i], v, offset, theta);
	}
	return;
}

/* magnitude()
 *  Calculated magnitude of vector *v.
 * arg *v: Vector
 * returns: Magnitude
 */
float magnitude(struct Vector *v) {
	float n;
	n = powf(v->x, 2.) + powf(v->y, 2.) + powf(v->z, 2.);
	n = sqrtf(n);
	return n;
}

/* normalise()
 *  Normalise vector *v in situ
 * arg *v: Vector to be normalised.
 * 
 * Error modes
 *  - If vector has 0 magnitude, will do nothing.
 */
void normalise(struct Vector *v) {
	float n = magnitude(v);
	if (n == 0)
		return;
	v->x /= n;
	v->y /= n;
	v->z /= n;
	return;
}

/* rotate_about()
 *  Will apply a rotation of theta radians about the bond *a-*b.
 *  Fixes atom *a's position, and rotates *b and all neighbours of *b.
 *  Should reset_check prior to running.
 * arg *a: Atom A (fixed!)
 * arg *b: Atom B
 * arg theta: angle in radians 
 */
int rotate_about(struct Atom *a, struct Atom *b, float theta) {
	// a is fixed.
	a->check = TRUE;
	a->rotate = FALSE;
	
	// offset is set to be a. So the system will be offset
	// such that a is at the origin.
	struct Vector offset;
	offset.x = a->v.x;
	offset.y = a->v.y;
	offset.z = a->v.z;
	
	struct Vector rotate_about;
	// we calculate rotate_about, which is (*b - *a)/|(*b - *a)|, e.g.
	// a normalised vector with a at the origin pointing to b.
	sub_vector(&(b->v), &(a->v), &rotate_about);
	normalise(&rotate_about);
	// and we apply rotation
	rotate(b, &rotate_about, &offset, theta);
	return 1;
}
	
/* bond_xyz()
 *  XYZ files do not contain bond information; this function adds bonds
 *  between all atoms where the distance between the atoms is < bl.
 * arg *m: Molecule pointer
 * arg bl: Maximum bond length
 * TODO: Make bl able to vary depending on atoms involved 
 *       e.g., C-H bond shorter than C-C.
 */
void bond_xyz(struct Molecule *m, float bl) {
	int i, j;
	struct Vector temp;
	for (i = 0; i < m->n_atoms; i++) {
		for (j = i+1; j < m->n_atoms; j++) {
			// calculate magnitude of vector joining them
			sub_vector(&(m->as[i].v), &(m->as[j].v), &temp);
			if (magnitude(&temp) < bl) {
				// and if below bl, add a bond.
				add_bond(&(m->as[i]), &(m->as[j]));
				printf("Bonding %d - %d\n", i, j);
			}
		}
	}
	return;
}

/* read_xyz()
 *  Reads xyz file *filename, and generates atomic positions and puts
 *  them into molecule *m. 
 * arg *m: Molecule pointer
 * arg *filename: Filename.
 * returns: -1 if file does not exist or is wrong format.
 *           1 if okay.
 */
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
			// first line in an xyz file gives the number of atoms.
			int k = sscanf(line, "%d\n", &(m->n_atoms));
			// if there is no number, close and exit.
			if (k != 1) {
				fclose(fp);
				return -1;
			} 
			// allocate memory for n_atoms atoms.
			m->as = (struct Atom *) malloc(sizeof(struct Atom) * m->n_atoms);
		} else if (c_line >= 2) {
			// c_line = 1 is a comment line, which we ignore.
			// if the current atom (c_line - 2) is greater than n_atoms - 1, 
			// then exit. eg if there are 3 atoms, then as[0-2].
			if (c_line - 2 > m->n_atoms - 1) {
				fclose(fp);
				return -1;
			}
			// ca is pointer to the current atom being operated on.
			ca = &(m->as[c_line - 2]);
			char *token = strtok(line, " \t");
			int i = 0;
			float x, y, z;
			char name[2];
			// while tokens are left...
			// i gives the column. First column is atom symbol, then x,y,z.
			while (token) {
				switch (i) {
					case 0: sprintf(name, "%2s", token); break;
					case 1: x = atof(token); break;
					case 2: y = atof(token); break;
					case 3: z = atof(token); break;
					// if there are more than 4 columns, crash and close.
					default: fclose(fp); return -1; break;
				}
				i++;

				token = strtok(NULL, " \t");
			}
			// add atom, and print out setup.
			add_atom(ca, x, y, z, name);
			printf("%s :: (%f, %f, %f) %d\n", ca->name, ca->v.x, ca->v.y, ca->v.z, c_line-2);
		}
		c_line++;
	}
	fclose(fp);
	return 1;
}

/* save_xyz()
 *  Outputs *m as an xyz file in *filename. Tab delimited, x y z are of 
 *  format +00.00000.
 * returns: -1 if file does not open
 *           1 otherwise
 */
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

/* print_dir()
 *  Prints molecule to screen
 */
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

/* run_script()
 *  Runs script in file *filename on molecule *m
 * Commands;
 *  - open FILENAME
 *     Opens xyz file at FILENAME
 *  - bond BL
 *     Generates bonds between atoms where distance < BL (angstroms)
 *  - rotate A B THETA
 *     Rotates atoms connected to B about the A-B axis an amount theta
 *     (radians)
 *  - output FILENAME
 *     Writes output as xyz file to FILENAME.
 * returns: -1 if file not available
 *           1 on success
 */
 
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
		// reset system checks
		reset_check(m);
        // read in line, up to point c
		sscanf(line, "%255s %n", command, &c);
		if (command[0] == '%')
			continue;
		if (strcmp(command, "open") == 0) {
			// if there is no second argument, complain.
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

/* main()
 *  Runs system. If a script is passed, runs the script. Else,
 *  enters a do-while loop operating on commands in an interactive manner
 * args SCRIPT: can pass script filename which will be run in.
 */
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
		/* Interactive mode commands
		 * - open FILENAME
		 *    opens xyz file FILENAME
		 * - bond BL
		 *    sets up bonds between all atoms with distance < BL angstroms
		 * - graph N
		 *    prints out graph starting at atom N
		 * - print
		 *    prints out entire system
		 * - rotate A B THETA
		 *    applies a rotation to atoms connected to B about the A-B axis
		 *    an amount theta radians.
		 * - output FILENAME
		 *    outputs as XYZ file into FILENAME
		 * - run SCRIPT
		 *    runs script in SCRIPT.
		 * - exit
		 *    exits program
		 */
		do {
			reset_check(&mol);
			printf("> ");
			scanf("%255s", command);
			printf("%s\n", command);
			if (strcmp(command, "open") == 0) {
				scanf("%255s", command);
				read_xyz(&mol, command);
			} else if (strcmp(command, "bond") == 0) {
				printf(" max bond length (ang) > ");
				scanf("%255s", command);
				bond_xyz(&mol, atof(command));
			} else if (strcmp(command, "graph") == 0) {
				printf(" start atom > ");
				if (scanf("%d", &n) == 1)
					print_moleculef(&(mol.as[n]), 0);
			} else if (strcmp(command, "print") == 0) {
				print_dir(&mol);
			} else if (strcmp(command, "rotate") == 0) {
				printf("This will perform a rotation about axis defined A -> B\n");
				printf("e.g. keeping A side atoms fixed, and rotating B connected atoms\n");
				printf("by an amount theta.\n");
				if (scanf("%d %d %f", &A, &B, &theta) != 3)
					continue;
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
	// clean up everything
	free_atoms(&mol);
	free(mol.as);
	return 1;
}
