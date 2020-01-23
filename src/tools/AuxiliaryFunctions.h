#pragma once

#include <string>
#include <vector>

#ifdef USE_INTEL_MKL
#include "mkl_pardiso.h"
#include "mkl_types.h"
#endif

#ifdef USE_PETSC
#include <petscksp.h>
#endif

class AuxiliaryFunctions
{
public:
	/***********************************************************************************************
	* \brief Constructor
	* \author Harikrishnan Sreekumar
	***********/
	AuxiliaryFunctions();
	/***********************************************************************************************
	* \brief Destructor
	* \author Harikrishnan Sreekumar
	***********/
	~AuxiliaryFunctions();
	/***********************************************************************************************
	* \brief Read a double vector from DAT file with format: <value>
	* \param[in] _fileName
	* \param[out] _vector::double
	* \author Harikrishnan Sreekumar
	***********/
	static std::vector<double> readDoubleVectorDatFormat(std::string _fileName);
    /***********************************************************************************************
	* \brief Read a int vector from DAT file with format: <value>
	* \param[in] _fileName
	* \param[out] _vector::double
	* \author Harikrishnan Sreekumar
	***********/
	static std::vector<int> readIntegerVectorDatFormat(std::string _fileName);
    /***********************************************************************************************
	* \brief Write a double vector to the DAT file with format: <value>
	* \param[in] _fileName
	* \param[in] _vector::double
	* \author Harikrishnan Sreekumar
	***********/
    static void writeDoubleVectorDatFormat(std::string _fileName, std::vector<double> &_vector);
#ifdef USE_INTEL_MKL
	/***********************************************************************************************
	* \brief Read a complex vector from DAT file with format: <value_real  value_imag>
	* \param[in] _fileName
	* \param[out] _vector::MKL_Complex16
	* \author Harikrishnan Sreekumar
	***********/
    static std::vector<MKL_Complex16> readComplexVectorDatFormat(std::string _fileName);
    /***********************************************************************************************
	* \brief [overloaded] Write a MKL_Complex16 vector to the DAT file with format: <value_real value_imag>
	* \param[in] _fileName
	* \param[in] _vector::MKL_Complex16
	* \author Harikrishnan Sreekumar
	***********/    
    static void writeComplexVectorDatFormat(std::string _fileName, std::vector<MKL_Complex16> &_vector);
    /***********************************************************************************************
	* \brief Converts a Complex Vector to Double Vector by ingoring the imaginary part
	* \param[in] MKL_Complex16 _vector
	* \param[out] _vector::double
	* \author Harikrishnan Sreekumar
	***********/
    static std::vector<double> convertComplexVectorToDouble(std::vector<MKL_Complex16> &_vector);
    /***********************************************************************************************
	* \brief Converts a Double Vector to Complex Vector zero as imaginary part
	* \param[in] double _vector
	* \param[out] _vector::MKL_Complex16
	* \author Harikrishnan Sreekumar
	***********/
    static std::vector<MKL_Complex16> convertDoubleVectorToComplex(std::vector<double> &_vector) ;
#endif
    
#ifdef USE_PETSC
    /***********************************************************************************************
	* \brief Read a complex vector from DAT file with format: <value_real  value_imag>
	* \param[in] _fileName
	* \param[out] _vector::PetscComplex
	* \author Harikrishnan Sreekumar
	***********/
    static std::vector<PetscComplex> readPetscComplexVectorDatFormat(std::string _fileName);
    /***********************************************************************************************
	* \brief [overloaded] Write a PetscComplex vector to the DAT file with format: <value_real value_imag>
	* \param[in] _fileName
	* \param[in] _vector::PetscComplex
	* \author Harikrishnan Sreekumar
	***********/    
    static void writeComplexVectorDatFormat(std::string _fileName, std::vector<PetscComplex> &_vector);
    
    static std::vector<PetscComplex> readPetscDoubleVectorDatFormatAsComplex(std::string _fileName);
#endif
};
