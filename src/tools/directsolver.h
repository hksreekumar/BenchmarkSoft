#pragma once

#ifdef USE_INTEL_MKL
#define MKL_DIRECT_CALL 1
#include "mkl.h"
#endif

#include <vector>


class directsolver{
public:
    /***********************************************************************************************
     * \brief Constructor
     * \author Harikrishnan Sreekumar
     ***********/
	directsolver();
    /***********************************************************************************************
     * \brief Constructor
     * \author Harikrishnan Sreekumar
     ***********/
#ifdef USE_INTEL_MKL
	directsolver(std::vector<MKL_INT> _ia, std::vector<MKL_INT> _ja, std::vector<double> _val, bool _sym);
#endif
    /***********************************************************************************************
     * \brief Constructor
     * \author Harikrishnan Sreekumar
     ***********/
#ifdef USE_INTEL_MKL
	directsolver(std::vector<MKL_INT> _ia, std::vector<MKL_INT> _ja, std::vector<MKL_Complex16> _val, bool _sym);
#endif
    /***********************************************************************************************
     * \brief Destructor
     *
     * \author Harikrishnan Sreekumar
     ***********/
	~directsolver();
    
    void initializeDirectSolver();
    
    void deinitializeDirectSolver();
    
    void factorize() ;
    
    void solve(std::vector<double> &_rhs, std::vector<double> &_sol);
#ifdef USE_INTEL_MKL
    void solve(std::vector<MKL_Complex16> &_rhs, std::vector<MKL_Complex16> &_sol);
#endif
    void setNumberOfCores(int _cores);
private:
    
    std::vector<double> my_mat_val;
#ifdef USE_INTEL_MKL
    std::vector<MKL_Complex16> my_mat_val_c;
    std::vector<MKL_INT> my_mat_ia;
    std::vector<MKL_INT> my_mat_ja;
#endif
    
    bool my_mat_sym;
    bool my_solver_complex;
#ifdef USE_INTEL_MKL
    // MKL Simple test
	/* Matrix data. */
	MKL_INT n;
	MKL_INT mtype;       /* Real symmetric matrix */
	MKL_INT nrhs;     /* Number of right hand sides. */
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;

	/* Auxiliary variables. */
	MKL_INT i;
	double ddum;          /* Double dummy */
	MKL_Complex16 ddum_c;
	MKL_INT idum;         /* Integer dummy. */
#endif
};
