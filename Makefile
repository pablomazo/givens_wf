FC=gfortran
FFLAGS=-O3

givens_wf: jacobi_diag.f90 givens_wf.f90 potential.f90
	$(FC) $(FFLAGS) -o givens.x $^
