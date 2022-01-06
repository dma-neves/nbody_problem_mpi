/**
 * Simulation over time of the N-Body problem
 * We don't pretend to follow Physics laws,
 * just some computation workload to play with.
 *
 * Based on the Reduced N-Body solver from:
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
vec *pos; /* just a 2D simulation */
vec *vel;

vec *loc_pos;
vec *loc_vel;
vec *loc_force;
vec *tmp_pos;
vec *tmp_force;
vec* tmp_pos_force;

int rank, size, locn;

int global_index(int p, int owner) { 
    
    return owner + p*size;
}

double get_mas(int owner, int p) {

    return mas[global_index(p, owner)];
}

int mod(int a, int b) {

    int r = a % b;
    return r < 0 ? r + b : r;
}

void initParticlesWithStride()
{
    pos = malloc(N_BODIES * sizeof(vec));
    vel = malloc(N_BODIES * sizeof(vec));

    for (int p = 0; p < N_BODIES; p++)
    {
        // Stride index that allows scattering/gathering the pos and vel arrays directly
        int stride_p = (p%size)*locn + p/size;

#ifdef DEBUG
        if(stride_pos >= N_BODIES)
            printf("ERROR un: Invalid index\n");
#endif

        mas[p] = 1;                /* 1 mass unit */
        pos[stride_p].x = drand48() * 100; /* 100x100 space */
        pos[stride_p].y = drand48() * 100;
        vel[stride_p].x = 0;
        vel[stride_p].y = 0;
    }
}

void printParticles(FILE *f)
{
    // print particles to text file f
    fprintf(f, "#  pos     vel\n");
    for (int p = 0; p < N_BODIES; p++)
    {
        // Get p's position and velocity from their stride location
        int stride_p = (p%size)*locn + p/size;

        fprintf(f, "%d (%g,%g) (%g,%g)\n",
                p, pos[stride_p].x, pos[stride_p].y, vel[stride_p].x, vel[stride_p].y);
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

void computeLocalForces(int q)
{
    for (int k = q+1; k < locn; k++)
    {
        vec force_qk = newton_sec_law(
            get_mas(rank,q), 
            loc_pos[q],
            get_mas(rank, k),
            loc_pos[k]
        );

        loc_force[q].x += force_qk.x;
        loc_force[q].y += force_qk.y;
        tmp_force[k].x -= force_qk.x;
        tmp_force[k].y -= force_qk.y;
    }
}

void computeReceivingForces(int owner, int q)
{
    int k = (owner <= rank) ? q+1 : q;
    for(; k < locn; k++) {

        vec force_qk = newton_sec_law(
            get_mas(rank,q), 
            loc_pos[q],
            get_mas(owner, k),
            tmp_pos[k]
        );

        loc_force[q].x += force_qk.x;
        loc_force[q].y += force_qk.y;
        tmp_force[k].x -= force_qk.x;
        tmp_force[k].y -= force_qk.y;
    }
}

void moveParticle(int q, double deltat)
{
#ifdef DEBUG
    if(loc_forcex[q] == 0 || loc_forcey[q] == 0)
        printf("Error: division by zero\n");
#endif
        
    loc_pos[q].x += deltat * loc_vel[q].x;
    loc_pos[q].y += deltat * loc_vel[q].y;
    loc_vel[q].x += deltat / get_mas(rank, q) * loc_force[q].x;
    loc_vel[q].y += deltat / get_mas(rank, q) * loc_force[q].y;
}

void simulateStep(double deltat)
{
    int source = mod(rank+1, size);
    int dest = mod(rank-1, size);

    memcpy(tmp_pos, loc_pos, locn*sizeof(vec));
    memset(loc_force, 0, locn*sizeof(vec));
    memset(tmp_force, 0, locn*sizeof(vec));

    for (int q = 0; q < locn; q++)
        computeLocalForces(q);

    for(int phase = 1; phase < size; phase++) {

        MPI_Sendrecv_replace(tmp_pos_force, 4*locn, MPI_DOUBLE, dest, rank, source, source, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int owner = (rank + phase) % size;

        for (int q = 0; q < locn; q++)
            computeReceivingForces(owner, q);
    }

    MPI_Sendrecv_replace(tmp_pos_force, 4*locn, MPI_DOUBLE, dest, rank, source, source, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for(int q = 0; q < locn; q++) {

        loc_force[q].x += tmp_force[q].x;
        loc_force[q].y += tmp_force[q].y;
    }

    for (int q = 0; q < locn; q++)
        moveParticle(q, deltat);
}

int main(int argc, char *argv[])
{

    /* ------------ args ------------ */

    int nSteps = 100; // default (you can give this at the command line)
    double time = 100;

    if (argc == 2)
    {
        nSteps = atoi(argv[1]); // number of steps
    }
    double deltat = time / nSteps;

    /* ------------ MPI init ------------ */

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == ROOT)
        printf("Started %d steps!\n", nSteps);

    /* ------------ create local arrays ------------ */

    locn = N_BODIES / size;
    loc_pos = malloc(locn * sizeof(vec));
    loc_vel = malloc(locn * sizeof(vec));
    loc_force = malloc(locn * sizeof(vec));
    tmp_pos_force = malloc(2 * locn * sizeof(vec));
    tmp_pos = tmp_pos_force;
    tmp_force = tmp_pos_force + locn;

    /* ------------ init -> scatter -> solve -> gather ------------ */

    if (rank == ROOT)
        initParticlesWithStride();

    clock_t t = clock();

    MPI_Bcast(mas, N_BODIES, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Scatter(pos, 2*locn, MPI_DOUBLE, loc_pos, 2*locn, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Scatter(vel, 2*locn, MPI_DOUBLE, loc_vel, 2*locn, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    for (int s = 0; s < nSteps; s++)
        simulateStep(deltat);

    MPI_Gather(loc_pos, 2*locn, MPI_DOUBLE, pos, 2*locn, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Gather(loc_vel, 2*locn, MPI_DOUBLE, vel, 2*locn, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    t = clock() - t;

     /* ------------ free local arrays ------------ */

    free(loc_pos);
    free(loc_vel);
    free(loc_force);
    free(tmp_pos_force);

    /* ------------ output result ------------ */

    if (rank == ROOT) {

        printf("time: %f s\n", t / (double)CLOCKS_PER_SEC);
        FILE* f = fopen("reduced_res", "w+");
        printParticles(f); // check if this solution is correct

        free(vel);
        free(pos);
    }

    MPI_Finalize();
    return 0;
}
