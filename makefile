MAKE=make
SRC=src
BIN=bin
DATA=data

all:
	$(MAKE) clean;
	$(MAKE) COMP_TYPE=STD -C $(SRC);
	mv $(SRC)/LibreGrowth $(BIN)/LibreGrowth;

clean:
	rm -f $(BIN)/LibreGrowth;
	rm -f $(DATA)/*;
	mkdir -p $(BIN);
	mkdir -p $(DATA);
	$(MAKE) clean -C $(SRC);

run: 
	$(MAKE) clean;
	$(MAKE) COMP_TYPE=STD -C $(SRC);
	mv $(SRC)/LibreGrowth $(BIN)/LibreGrowth;
	ulimit -s unlimited; cd ./$(DATA); time ../$(BIN)/LibreGrowth > log.dat;

run-omp: 
	$(MAKE) clean;
	$(MAKE) COMP_TYPE=OMP -C $(SRC);
	mv $(SRC)/LibreGrowth $(BIN)/LibreGrowth;
	export OMP_NUM_THREADS=8; ulimit -s unlimited; cd ./$(DATA); time ../$(BIN)/LibreGrowth > log.dat;

run-dbg:
	$(MAKE) clean;
	$(MAKE) COMP_TYPE=DBG -C $(SRC);
	mv $(SRC)/LibreGrowth $(BIN)/LibreGrowth;
	cd ./$(DATA); ddd ../$(BIN)/LibreGrowth;


