/**
 * Simulation over time of the N-Body problem
 * We don't pretend to follow Physics laws,
 * just some computation workload to play with.
 *
 * Based on the Basic N-Body solver from:
 * An introduction to parallel programming, 2nd Edition
 * Peter Pacheco, Matthew Malensek
 *
 * Students: David Neves 55539 and Rodrigo Mesquita 55902
 * Teacher: Vitor Duarte FCT/UNL 2021
 * CAD - 2021/2022
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>

//#define N_BODIES    2048
#define N_BODIES 64
#define Gconst 0.1
#define ROOT 0

typedef struct _vec {

    double x;
    double y;
} vec;

double mas[N_BODIES];
vec pos[N_BODIES]; /* just a 2D simulation */

vec *vel;
vec *loc_pos;
vec *loc_vel;
vec *loc_force;

int rank, size, locn;

void initParticles()
{

    vel = malloc(N_BODIES * sizeof(vec));

    for (int p = 0; p < N_BODIES; p++)
    {
        mas[p] = 1;                /* 1 mass unit */
        pos[p].x = drand48() * 100; /* 100x100 space */
        pos[p].y = drand48() * 100;
        vel[p].x = 0;
        vel[p].y = 0;
    }
}

void printParticles(FILE *f)
{
    // print particles to text file f
    fprintf(f, "#  pos     vel\n");
    for (int p = 0; p < N_BODIES; p++)
    {
        fprintf(f, "%d (%g,%g) (%g,%g)\n",
                p, pos[p].x, pos[p].y, vel[p].x, vel[p].y);
    }
}

vec newton_sec_law(double mas_a, vec pos_a, double mas_b, vec pos_b) {

    double xdiff = pos_a.x - pos_b.x;
    double ydiff = pos_a.y - pos_b.y;

    double dist = sqrt(xdiff * xdiff + ydiff * ydiff);
    double distCub = dist * dist * dist;

    return (vec) {

        -Gconst * mas_a * mas_b / distCub * xdiff,
        -Gconst * mas_a * mas_b / distCub * ydiff
    };
}

void computeForces(int q, int global_q)
{
    // based on the basic solver from book (not the reduced solver)
    for (int k = 0; k < N_BODIES; k++)
    {
        if (k == global_q)
            continue; // ignore itself

        vec force_qk = newton_sec_law(

            mas[q],
            loc_pos[q],
            mas[k],
            pos[k]
        );

        loc_force[q].x += force_qk.x;
        loc_force[q].y += force_qk.y;
    }
}

void moveParticle(int q, double deltat)
{
    loc_pos[q].x += deltat * loc_vel[q].x;
    loc_pos[q].y += deltat * loc_vel[q].y;
    loc_vel[q].x += deltat / mas[q] * loc_force[q].x;
    loc_vel[q].y += deltat / mas[q] * loc_force[q].y;
}

void simulateStep(double deltat)
{

    memset(loc_force, 0, sizeof(vec) * locn);

    for (int q = 0; q < locn; q++)
        computeForces(q, rank*locn+q);

    for (int q = 0; q < locn; q++)
        moveParticle(q, deltat);

    MPI_Allgather(loc_pos, locn*2, MPI_DOUBLE, pos, locn*2, MPI_DOUBLE, MPI_COMM_WORLD);
}

int main(int argc, char *argv[])
{

    int nSteps = 100; // default (you can give this at the command line)
    double time = 100;

    if (argc == 2)
    {
        nSteps = atoi(argv[1]); // number of steps
    }
    double deltat = time / nSteps;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == ROOT)
        printf("Started %d steps!\n", nSteps);

    locn = N_BODIES / size;
    loc_pos = pos + rank*locn;
    loc_vel = malloc(locn * sizeof(vec));
    loc_force = malloc(locn * sizeof(vec));

    if (rank == ROOT)
        initParticles();

    clock_t t = clock();

    MPI_Bcast(mas, N_BODIES, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(pos, 2*N_BODIES, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Scatter(vel, 2*locn, MPI_DOUBLE, loc_vel, 2*locn, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    for (int s = 0; s < nSteps; s++)
        simulateStep(deltat);

    MPI_Gather(loc_vel, 2*locn, MPI_DOUBLE, vel, 2*locn, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    t = clock() - t;

    if (rank == ROOT)
    {
        printf("time: %f s\n", t / (double)CLOCKS_PER_SEC);

        FILE* f = fopen("basic_res", "w+");
        printParticles(f);
    }

    MPI_Finalize();
    return 0;
}
