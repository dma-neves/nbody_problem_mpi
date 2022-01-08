LIBS=-lcsfml-graphics -lcsfml-window -lcsfml-system -lcsfml-audio -lm

all: nbody_basic nbody_reduced nbody_serial nbody_renderer

nbody_serial:	nbody_serial.c
	cc -o $@ $< -lm

nbody_reduced: nbody_reduced.c
	mpicc -o $@ $< -lm

nbody_basic: nbody_basic.c
	mpicc -o $@ $< -lm

nbody_renderer: nbody_renderer.c
	gcc -o $@ $< $(LIBS)

runreduced:
	mpirun -n 4 nbody_reduced

runbasic:
	mpirun -n 4 nbody_basic

testreduced:
	./nbody_serial
	mpirun -n 4 nbody_reduced
	diff serial_res reduced_res

testbasic:
	./nbody_serial
	mpirun -n 4 nbody_basic
	diff serial_res basic_res

clean:
	rm -f nbody_serial nbody_reduced nbody_basic nbody_renderer
	rm serial_res basic_res reduced_res
