#include <iostream>
#include <stdio.h>

#ifdef USE_INTEL_MKL
#define MKL_DIRECT_CALL 1
#include "mkl.h"
#endif

#ifdef USE_INTEL_MKL
#include "mkl_pardiso.h"
#include "mkl_types.h"
#endif

#ifdef USE_PETSC
#include <petscksp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <MathTools.h>
#include <AuxiliaryFunctions.h>
#include <chrono>

#include <directsolver.h>

#include <sstream>

int main(int argc,char **argv){
    std::cout<<"Starting to Benchmark..."<<std::endl;
    
    std::cout<<std::endl<<"---- Environment -------------------------------------------------------------------------"<< std::endl;
    #ifdef _CXX_VER
        #define CXX_VER _CXX_VER
    #endif
    #ifdef _CXX_ID
        #define CXX_ID _CXX_ID
    #endif
    std::cout << "  | Compiler     : " << CXX_ID << " version "<< CXX_VER << "." << std::endl;
    int len = 198;
    char buf[198];
#ifdef USE_PETSC
    PetscGetVersion(buf, len);
    std::cout << "  | PETSC        : " << buf << std::endl;
#else
    std::cout << "  | PETSC        : [NOT DETECTED]" << std::endl;
#endif
#ifdef USE_MPI
    int _mp = MPI_Get_library_version(buf, &len);
    std::cout << "  | MPI          : " << buf << std::endl;
#else
    std::cout << "  | MPI          : [NOT DETECTED]" << std::endl;    
#endif
#ifdef USE_INTEL_MKL
    // Print MKL Version
    mkl_get_version_string(buf, len);
    std::cout << "  | MKL          : " << buf << std::endl;
#else
    std::cout << "  | MKL          : [NOT DETECTED]" << std::endl;
#endif
        
    std::cout<<"------------------------------------------------------------------------------------------"<< std::endl << std::endl;
        
	#ifdef USE_INTEL_MKL
    /*---------------------------*/
    bool IS_SYMMETRIC = false;
    bool IS_COMPLEX = true;
    bool IS_HACK = false;
    int  NUM_THREADS = 6;
    /*----------------------------*/
    
    std::vector<double> _mat_val;
    std::vector<MKL_Complex16> _mat_val_c;
    std::vector<int> _mat_ia_D;
    std::vector<int> _mat_ja_D;
    std::vector<double> _rhs;
    std::vector<double> _sol;
    std::vector<MKL_Complex16> _rhs_c;
    std::vector<MKL_Complex16> _sol_c;
    std::string iopath = "../iodir/";
    /* Case: Jo */
    //std::string matprefix = "Test1e5_KDYN_UNSY";
    //std::string matprefix = "KK_n_jo";
    //std::string rhsname = matprefix + "_rhs.dat"; 
    
    /* Case: PlateSuperFine 261004 */
    std::string matprefix = "m_K";
    std::string rhsname = matprefix + "_rhs.dat"; 
    
    /* Case: PlateSuperFine visc_m */
    //std::string matprefix = "visc_m_K";
    //std::string rhsname = matprefix + "_rhs.dat"; 
    
    
    _mat_ia_D  = AuxiliaryFunctions::readIntegerVectorDatFormat(iopath + "input/" + matprefix + "_csr_ia.dat");
    _mat_ja_D  = AuxiliaryFunctions::readIntegerVectorDatFormat(iopath + "input/" + matprefix + "_csr_ja.dat");
    // Type cast to MKL_INT
    std::vector<MKL_INT> _mat_ia(_mat_ia_D.begin(), _mat_ia_D.end());
    std::vector<MKL_INT> _mat_ja(_mat_ja_D.begin(), _mat_ja_D.end());
    if(!IS_HACK) {
        if(IS_COMPLEX){
            _mat_val_c = AuxiliaryFunctions::readComplexVectorDatFormat(iopath + "input/" + matprefix + "_csr_val.dat");
            _rhs_c     = AuxiliaryFunctions::readComplexVectorDatFormat(iopath + "input/" + rhsname);
            _sol_c.resize(_mat_ia.size()-1);
        }
        else {
            _mat_val = AuxiliaryFunctions::readDoubleVectorDatFormat(iopath + "input/" + matprefix + "_csr_val.dat");
            _rhs     = AuxiliaryFunctions::readDoubleVectorDatFormat(iopath + "input/" + rhsname);
            _sol.resize(_mat_ia.size()-1);
        }
    } else {
        if(IS_COMPLEX){
            _mat_val_c = AuxiliaryFunctions::readComplexVectorDatFormat(iopath + "input/" + matprefix + "_csr_val.dat");
            _rhs_c     = AuxiliaryFunctions::readComplexVectorDatFormat(iopath + "input/" + rhsname);
            _sol.resize(_mat_ia.size()-1);
            
            _mat_val = AuxiliaryFunctions::convertComplexVectorToDouble(_mat_val_c);
            _rhs = AuxiliaryFunctions::convertComplexVectorToDouble(_rhs_c);
            
            _mat_val_c.clear();
            _rhs_c.clear();
            IS_COMPLEX = false;
        }
    }
    
    
    _mat_ia_D.clear(); _mat_ja_D.clear();
    
    
    std::cout<<"Input read: ia "<< _mat_ia.size() << " , ja: " << _mat_ja.size() << " , val: "<< _mat_val.size() << std::endl;
    
    directsolver* pardisoSolver;
    if(IS_COMPLEX)
        pardisoSolver = new directsolver(_mat_ia, _mat_ja, _mat_val_c, IS_SYMMETRIC);
    else
        pardisoSolver = new directsolver(_mat_ia, _mat_ja, _mat_val, IS_SYMMETRIC);
        
    pardisoSolver->setNumberOfCores(NUM_THREADS);
    pardisoSolver->initializeDirectSolver();
    
    std::chrono::high_resolution_clock::time_point startTime;
    std::chrono::high_resolution_clock::time_point stopTime;
    startTime = std::chrono::high_resolution_clock::now();
    pardisoSolver->factorize();
    if(IS_COMPLEX)
        pardisoSolver->solve(_rhs_c,_sol_c);
    else
        pardisoSolver->solve(_rhs,_sol);
    
	printf("\nSolve completed! ");
    stopTime = std::chrono::high_resolution_clock::now();
    std::cout<<"\n\nPARDISO: Performance | " << std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime).count() << " micro sec | " << std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count() << " milli sec | " << std::chrono::duration_cast<std::chrono::seconds>(stopTime - startTime).count() << " sec | " << std::endl;
    
    pardisoSolver->deinitializeDirectSolver();
	printf("\nThe solution of the system is: ");
	for (int i = 0; i < 5; i++)
	{
        if(IS_COMPLEX)
            std::cout<<"\n x["<<i<<"] = "<< _sol_c[i].real << " + i "<< _sol_c[i].imag;
        else
            std::cout<<"\n x["<<i<<"] = "<< _sol[i];
	}
	printf("\n");
    
    if(IS_COMPLEX)
        AuxiliaryFunctions::writeComplexVectorDatFormat(iopath + "output/" + matprefix + "_SOL_MKL.dat", _sol_c);
    else
        AuxiliaryFunctions::writeDoubleVectorDatFormat(iopath + "output/" + matprefix + "_SOL_MKL.dat", _sol);
#endif
    
#ifdef USE_PETSC
    
    /*---------------------------*/
    bool PETSC_IS_SYMMETRIC = false;
    bool PETSC_IS_COMPLEX = true;
    bool PETSC_IS_HACK = false;
    int  PETSC_NUM_THREADS = 1;
    /*----------------------------*/
    
    std::vector<PetscComplex> _p_mat_val_c;
    std::vector<int> _p_mat_ia_D;
    std::vector<int> _p_mat_ja_D;
    
    std::vector<PetscComplex> _p_rhs_c;
    std::vector<PetscComplex> _p_sol_c;
    
    std::string p_iopath = "../iodir/";
    /* Case: Jo */
    //std::string matprefix = "Test1e5_KDYN_UNSY";
    //std::string matprefix = "KK_n_jo";
    //std::string rhsname = matprefix + "_rhs.dat"; 
    
    /* Case: PlateSuperFine 261004 */
    std::string p_matprefix = "m_K";
    std::string p_rhsname = p_matprefix + "_rhs.dat"; 
    
    _p_mat_ia_D  = AuxiliaryFunctions::readIntegerVectorDatFormat(p_iopath + "input/" + p_matprefix + "_csr_ia.dat");
    
    _p_mat_ja_D  = AuxiliaryFunctions::readIntegerVectorDatFormat(p_iopath + "input/" + p_matprefix + "_csr_ja.dat");
    // Type cast to PetscInt
    std::vector<PetscInt> _p_mat_ia(_p_mat_ia_D.begin(), _p_mat_ia_D.end());
    std::vector<PetscInt> _p_mat_ja(_p_mat_ja_D.begin(), _p_mat_ja_D.end());
    if(!PETSC_IS_HACK) {
        if(PETSC_IS_COMPLEX){
            _p_mat_val_c = AuxiliaryFunctions::readPetscComplexVectorDatFormat(p_iopath + "input/" + p_matprefix + "_csr_val.dat");
            _p_rhs_c     = AuxiliaryFunctions::readPetscComplexVectorDatFormat(p_iopath + "input/" + p_rhsname);
            //_p_sol_c.resize(_p_mat_ia.size()-1);
        }
        else {
            exit(EXIT_FAILURE);
        }
    } else {
        if(PETSC_IS_COMPLEX){
            exit(EXIT_FAILURE);
        }
    }
    
    _p_mat_ia_D.clear(); _p_mat_ja_D.clear();
    
    std::cout<<"Input read: ia "<< _p_mat_ia.size() << " , ja: " << _p_mat_ja.size() << " , val: "<< _p_mat_val_c.size() << std::endl;
    PetscInt m_NumberOfUnknowns = _p_mat_ia.size() -1 ;
    
    /* -- Convert to zero based -- */
    for(int i_ia = 0; i_ia < _p_mat_ia.size() ; i_ia++) {
        _p_mat_ia[i_ia] = _p_mat_ia[i_ia] - 1;
    }
    for(int i_ja = 0; i_ja < _p_mat_ja.size() ; i_ja++) {
        _p_mat_ja[i_ja] = _p_mat_ja[i_ja] - 1;
    }
    
    PetscErrorCode p_error;
    Mat K_dy, F;
    PetscBool      nonzeroguess = PETSC_FALSE,changepcside = PETSC_FALSE;
    PetscInt       i,col[3],its;
    PetscMPIInt    size;
    char help[] = "Solves a tridiagonal linear system with KSP.\n\n";
    
    p_error = PetscInitialize(&argc,&argv,(char*)0,help);if (p_error) return p_error;
    p_error = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(p_error);
    //if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
    p_error = PetscOptionsGetInt(NULL,NULL,"-n",&m_NumberOfUnknowns,NULL);CHKERRQ(p_error);
    p_error = PetscOptionsGetBool(NULL,NULL,"-nonzero_guess",&nonzeroguess,NULL);CHKERRQ(p_error);
  
    //p_error = MatCreateMPIAIJWithArrays( PETSC_COMM_WORLD , m_NumberOfUnknowns, m_NumberOfUnknowns,PETSC_DECIDE, PETSC_DECIDE, &_p_mat_ia[0], &_p_mat_ja[0], &_p_mat_val_c[0], &K_dy );
    p_error = MatCreateSeqAIJWithArrays( PETSC_COMM_WORLD , m_NumberOfUnknowns, m_NumberOfUnknowns, &_p_mat_ia[0], &_p_mat_ja[0], &_p_mat_val_c[0], &K_dy );
    
    its = 0;

    std::cout<< " > Solving linear system of equations using PETSC... " << std::endl;
    //resetPreconditioner = true;

    KSP  m_KSP;         // PETSc's linear equation solver
    PC   m_PC;          // Prrconditioner
    Vec  m_x;           // Solution
    Vec  m_RHS;         // RHS
    
    MatFactorInfo  info;
    IS             perm,iperm;
    
    p_error = KSPCreate(PETSC_COMM_WORLD, &m_KSP); 
    p_error = KSPGetPC(m_KSP, &m_PC); 

    PetscInt solverSetup = 1;
    PetscBool found = PETSC_FALSE;
    PetscOptionsGetInt(PETSC_NULL,"", "-solver", &solverSetup, &found);
    
    switch(solverSetup){
        case 1:
            std::cout<< " > Initializing direct mumps solver "<<std::endl;
            p_error = KSPSetType(m_KSP, KSPPREONLY); 
            p_error = PCSetType(m_PC, PCLU); 
            p_error = PCFactorSetMatSolverType(m_PC,"mumps"); 
            break;
        case 2:
            std::cout<< " > Initializing direct pardiso solver "<<std::endl;
            p_error = KSPSetType(m_KSP, KSPPREONLY); 
            p_error = PCSetType(m_PC, PCLU); 
            p_error = PCFactorSetMatSolverType(m_PC,"mkl_pardiso"); 
            break;

            break;
        case 4:
            std::cout<<"  initializing direct mumps solver (default / case 4) \n";
            p_error = KSPSetType(m_KSP, KSPPREONLY); 
            p_error = PCSetType(m_PC, PCLU); 
            p_error = PCFactorSetMatSolverType(m_PC,"mkl_pardiso"); 
            solverSetup = 1;
            break;
        default:
            std::cerr<<" > Invalid Solver Type!"<< std::endl;
            exit(EXIT_FAILURE);
    }
    
    p_error = VecCreate(PETSC_COMM_WORLD, &m_x); 
    p_error = VecSetSizes(m_x, PETSC_DECIDE, m_NumberOfUnknowns); 

    p_error = VecSetFromOptions(m_x); 
    p_error = PetscObjectSetName((PetscObject)m_x ,"x"); 
    
    VecDuplicate(m_x, &m_RHS);
    VecPlaceArray(m_RHS, &_p_rhs_c[0]);
    /*
    PetscInt size_x, size_F;
    VecGetSize(m_x,&size_x);
    VecGetSize(m_RHS,&size_F);
    
    std::cout<<"Size of x "<< size_x << " F: "<< size_F << std::endl;*/

    p_error = PetscObjectSetName((PetscObject)m_RHS ,"F"); 
    
    // This loop is placed here AND in funtion "initializeSolverObjects" as some setup require KSPSetOperators before (e.g. mumps icntl)
    std::chrono::high_resolution_clock::time_point p_startTime;
    std::chrono::high_resolution_clock::time_point p_stopTime;
    p_startTime = std::chrono::high_resolution_clock::now();
    
    p_error = KSPSetOperators(m_KSP, K_dy, K_dy); 
    p_error = KSPSolve(m_KSP, m_RHS, m_x); 
    p_error = KSPGetIterationNumber(m_KSP, &its); 

    printf("\nSolve completed! ");
    p_stopTime = std::chrono::high_resolution_clock::now();
    std::cout<<"\n\nPetSc  : Performance | " << std::chrono::duration_cast<std::chrono::microseconds>(p_stopTime - p_startTime).count() << " micro sec | " << std::chrono::duration_cast<std::chrono::milliseconds>(p_stopTime - p_startTime).count() << " milli sec | " << std::chrono::duration_cast<std::chrono::seconds>(p_stopTime - p_startTime).count() << " sec | " << std::endl;
    
    printf("\nThe solution of the system is: ");
    
    PetscComplex* _array;
    if(PETSC_COMPLEX){
        p_error = VecGetArray(m_x,&_array);
        _p_sol_c.insert(_p_sol_c.end(), &_array[0], &_array[m_NumberOfUnknowns]);
    }
    else
        exit(EXIT_FAILURE);
        
    for (int i = 0; i < 5; i++)
    {
        if(PETSC_IS_COMPLEX){
            std::cout<<"\n x["<< i <<"] = "<< _p_sol_c[i].real() << " + i "<< _p_sol_c[i].imag();
        }
        else
            exit(EXIT_FAILURE);
    }
    printf("\n");
    
    p_error = VecDestroy(&m_x); 
    p_error = VecDestroy(&m_RHS);
    p_error = MatDestroy(&K_dy);   
    p_error = KSPDestroy(&m_KSP);
    
    if(PETSC_IS_COMPLEX)
        AuxiliaryFunctions::writeComplexVectorDatFormat(p_iopath + "output/" + p_matprefix + "_SOL_PETSC.dat", _p_sol_c);
    else
        exit(EXIT_FAILURE);
#endif
    return 0;
}
