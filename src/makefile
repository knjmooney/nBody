CC        = mpicc
PROJ      = n-body
OBJECTS   = $(PROJ).o point.o integrator.o grid.o para_grid.o
CFLAGS    = -g -O2 -fdiagnostics-color=auto -Wall -Wextra -pedantic -std=gnu99
LDLIBS    = -lm

# Arguments for making graphs and animations
DATAFILE = data/animation.dat
TESTARGS = -s 10000000 -t 0.000001 -n 6 -w 4 -h 4

PARA_DATAFILE = data_par/parafile_0.dat 

$(PROJ): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDLIBS)

test: $(PROJ)
	 ./$(PROJ) $(TESTARGS)

.PHONY: clean		
clean:
	rm -f $(PROJ) $(OBJECTS)

$(DATAFILE): $(PROJ)
	./$(PROJ) $(TESTARGS) > $(DATAFILE)

$(PARA_DATAFILE): $(PROJ)
	# mpirun -n 9 ./$(PROJ) -n 18 -s 1000 -w 5 -h 5 -t 0.0001
	sbatch slurm.sh

animation: $(DATAFILE)
	cd data; gnuplot animation.gp

para_animation: $(PARA_DATAFILE)
	cd data_par; gnuplot plot.gp


energy: $(DATAFILE)
	cd data; gnuplot energy.gp
