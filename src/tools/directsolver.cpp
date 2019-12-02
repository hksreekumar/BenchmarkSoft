#include <directsolver.h>
#include <iostream>

directsolver::directsolver(){
}

#ifdef USE_INTEL_MKL
directsolver::directsolver(std::vector<MKL_INT> _ia, std::vector<MKL_INT> _ja, std::vector<double> _val, bool _sym){
    my_mat_ia = _ia;
    my_mat_ja = _ja;
    my_mat_val = _val;
    
    my_mat_sym = _sym;
    
    n = _ia.size()-1;
    if(_sym) {
        mtype = -2;
        std::cout<<">> Solver set to real symmetric"<<std::endl;
    } else {
        mtype = 11;
        std::cout<<">> Solver set to real unsymmetric"<< std::endl;
    }
    nrhs = 1;
    my_solver_complex = false;
    std::cout<< ">> MKL running on #" << mkl_get_max_threads() << " OpenMPI cores. Change if necessary."<< std::endl;
}
#endif

#ifdef USE_INTEL_MKL
directsolver::directsolver(std::vector<MKL_INT> _ia, std::vector<MKL_INT> _ja, std::vector<MKL_Complex16> _val, bool _sym){
    my_mat_ia = _ia;
    my_mat_ja = _ja;
    my_mat_val_c = _val;
    
    my_mat_sym = _sym;
    
    n = _ia.size()-1;
    if(_sym) {
        mtype = 6;
        std::cout<<">> Solver set to complex symmetric"<<std::endl;
    } else {
        mtype = 13;
        std::cout<<">> Solver set to complex unsymmetric"<< std::endl;
    }
    nrhs = 1;
    my_solver_complex = true;
}
#endif

directsolver::~directsolver() {
}

void directsolver::initializeDirectSolver(){
#ifdef USE_INTEL_MKL
    /* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
	pardisoinit(pt, &mtype, iparm);
    /* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
	{
		//iparm[i] = 0;
	}
	
	iparm[0] = 1;    // No solver defaults
    iparm[1] = 2;    // Fill-in reordering from METIS 
    iparm[7] = 2;    // Max numbers of iterative refinement steps
	iparm[9] = 13;   // Perturb the pivot elements with 1E-13
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	iparm[23] = 1;   // 2-level factorization
	iparm[35] = 0; 
	iparm[36] = -99; // VBSR format
	iparm[56] = 2;
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[17] = -1;	 // Output: Number of nonzeros in the factor LU
	iparm[18] = -1;	 // Output: Report Mflops
	iparm[19] = 0;	 // Output: Number of CG iterations
	// iparm[27] = 1;   // PARDISO checks integer arrays ia and ja. In particular, PARDISO checks whether column indices are sorted in increasing order within each row.
	maxfct = 1;	    // max number of factorizations
	mnum = 1;		// which factorization to use
	msglvl = 1;		// do NOT print statistical information
	error = 1;		// Initialize error flag 
    nrhs = 1; // number of right hand side
#endif
}

void directsolver::factorize() {
#ifdef USE_INTEL_MKL
    /* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	if(my_solver_complex)
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &my_mat_val_c[0], &my_mat_ia[0], &my_mat_ja[0], &idum, &nrhs, iparm, &msglvl, &ddum_c, &ddum_c, &error);
    else
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &my_mat_val[0], &my_mat_ia[0], &my_mat_ja[0], &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %lld", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %lld", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %lld \n", iparm[18]);
	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	if(my_solver_complex)
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &my_mat_val_c[0], &my_mat_ia[0], &my_mat_ja[0], &idum, &nrhs, iparm, &msglvl, &ddum_c, &ddum_c, &error);
    else
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &my_mat_val[0], &my_mat_ia[0], &my_mat_ja[0], &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
        
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %lld", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");
#endif
}

void directsolver::solve(std::vector<double> &_rhs, std::vector<double> &_sol){
#ifdef USE_INTEL_MKL
    /* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
	
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,&n, &my_mat_val[0], &my_mat_ia[0], &my_mat_ja[0], &idum, &nrhs, iparm, &msglvl, &_rhs[0], &_sol[0], &error);
    
	if (error != 0)
	{
		printf("\nERROR during solution: %lld", error);
		exit(3);
	}
#endif
}
#ifdef USE_INTEL_MKL
void directsolver::solve(std::vector<MKL_Complex16> &_rhs, std::vector<MKL_Complex16> &_sol){
    /* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, &my_mat_val_c[0], &my_mat_ia[0], &my_mat_ja[0], &idum, &nrhs, iparm, &msglvl, &_rhs[0], &_sol[0], &error);
	if (error != 0)
	{
		printf("\nERROR during solution: %lld", error);
		exit(3);
	}
}
#endif

void directsolver::setNumberOfCores(int _cores) {
#ifdef USE_INTEL_MKL
    mkl_set_num_threads(_cores);
    std::cout<<">> Intel MKL set to #" << mkl_get_max_threads() <<" cores." << std::endl;
#endif
}

void directsolver::deinitializeDirectSolver() {
#ifdef USE_INTEL_MKL
    /* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1;           /* Release internal memory. */
	if(my_solver_complex)
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, &my_mat_ia[0], &my_mat_ja[0], &idum, &nrhs,iparm, &msglvl, &ddum_c, &ddum_c, &error);
    else
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, &my_mat_ia[0], &my_mat_ja[0], &idum, &nrhs,iparm, &msglvl, &ddum, &ddum, &error);
#endif
}
