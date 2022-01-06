#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <SFML/Graphics.h>

#define N_BODIES 16
#define Gconst  0.1

#define WIDTH 500
#define HEIGHT 500
#define WWIDTH 100
#define WHEIGHT 100
#define DT_MULT 20

int color[N_BODIES][3];
double mas[N_BODIES];
double posx[N_BODIES];  /* just a 2D simulation */
double posy[N_BODIES];
double velx[N_BODIES];
double vely[N_BODIES];
double forcex[N_BODIES];
double forcey[N_BODIES];

char dead[N_BODIES];

void init_sfml();
void render();
void handleEvents();

sfRenderWindow* window;
sfVideoMode mode;
sfEvent event;
int running = 1;
sfClock* sfclock;

void initParticles() {

    for (int p=0; p<N_BODIES; p++) {

        posx[p]=drand48()*WWIDTH/2 + WWIDTH/4;
        posy[p]=drand48()*WHEIGHT/2 + WHEIGHT/4;
        mas[p]=drand48(); if(mas[p] < 0.25) mas[p] = 0.25;
        velx[p]=(drand48()-0.5)*0.2;
        vely[p]=(drand48()-0.5)*0.2;
        color[p][0]=drand48()*255;
        color[p][1]=drand48()*255;
        color[p][2]=drand48()*255;
        dead[p] = 0;
    }
}

void computeForces(int q) {

    if(dead[q]) return;

    for (int k=q+1; k<N_BODIES; k++) {

        if(dead[k]) continue;
        
        double xdiff=posx[q] - posx[k];
        double ydiff=posy[q] - posy[k];
        double dist=sqrt(xdiff*xdiff+ydiff*ydiff);

        if(dist < 1) {

            // kill/consume smaller particle
            int i = mas[q] > mas[k] ? k : q;
            dead[i] = 1;
        }

        double distCub=dist*dist*dist;

        double force_qk_x = -Gconst * mas[q] * mas[k] / distCub * xdiff;
        double force_qk_y = -Gconst * mas[q] * mas[k] / distCub * ydiff;

        forcex[q] += force_qk_x;
        forcey[q] += force_qk_y;
        forcex[k] -= force_qk_x;
        forcey[k] -= force_qk_y;
    }
}

void moveParticle(int q, double deltat) {

    if(!dead[q]) {

        posx[q] += deltat*velx[q];
        posy[q] += deltat*vely[q];
        velx[q] += deltat/mas[q] * forcex[q];
        vely[q] += deltat/mas[q] * forcey[q];
    }
}

void simulateStep(double deltat) {

    memset(forcex, 0, sizeof forcex);
    memset(forcey, 0, sizeof forcey);
    for ( int q=0; q<N_BODIES; q++)
        computeForces(q);
    for (int q=0; q<N_BODIES; q++)
        moveParticle(q, deltat);
}

int main(int argc, char *argv[]) {

    init_sfml();
    initParticles();

    while(running)
    {
        float dt = sfTime_asSeconds( sfClock_restart(sfclock) )*DT_MULT;

        handleEvents();
        simulateStep(dt);
        render();
    }

    return 0;
}

/* ------------- sfml ------------- */

void init_sfml() {

    sfclock = sfClock_create();
    mode.height = 500;
    mode.width = 500;
    mode.bitsPerPixel = 32;
    window = sfRenderWindow_create(mode, "n-body", sfResize | sfClose, NULL);

}

void handleEvents()
{
    while(sfRenderWindow_pollEvent(window, &event))
    {
        if(event.type == sfEvtClosed || sfKeyboard_isKeyPressed(sfKeyEscape))
        {
            running = 0;
            sfRenderWindow_close(window);
        }
    }
}


void render() {

    //sfRenderWindow_clear(window, sfBlack);

    sfCircleShape* circle = sfCircleShape_create();

    for(int p = 0; p < N_BODIES; p++) {

        if(!dead[p]) {

            sfCircleShape_setFillColor(circle, sfColor_fromRGB(color[p][0],color[p][1],color[p][2]));
            // float rad = 5.0f*mas[p];
            // if(rad < 2) rad = 2;
            float rad = 2;
            sfCircleShape_setRadius(circle, rad);
            sfVector2f pos = { .x = posx[p]*WIDTH/WWIDTH - rad/2, .y = posy[p]*HEIGHT/WHEIGHT - rad/2 };
            sfCircleShape_setPosition(circle, pos);
            sfRenderWindow_drawCircleShape(window, circle, NULL);
        }
    }

    sfRenderWindow_display(window);
}