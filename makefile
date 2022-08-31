ALL: run_main 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

run_main: params_mod.o utilities_mod.o indexing_mod.o grid_mod.o soln_mod.o Au_mod.o Av_mod.o Ap_mod.o brinkman_mod.o Bu_mod.o Bv_mod.o initial_P_mod.o uRHS_mod.o vRHS_mod.o pRHS_mod.o nonlin_mod.o reshape_mod.o correction_mod.o fluids_mod.o main.o 
	-${FLINKER} -o $@ $^ ${PETSC_LIB}

CLEAN:
	 rm run_main *.mod *.o
