
#ifndef INFAM_OUTPUTH5_H
#define INFAM_OUTPUTH5_H

//#include "outputfile.h"
#include <vector>
#include "fileh5.h"

#ifdef HAVE_PETSC
#include "petscconf.h"
#include "petscksp.h"
#include "petscao.h"
#include "slepceps.h"
#endif // HAVE_PETSC

/**
* @brief Singleton H5 class supporting H5 write
* @author Harikrishnan Sreekumar
* @date 18.03.2020
*
* cOutputH5 writes data to a h5 file.
*/

class cOutputH5: public cFileH5
{
public:
    /**********************************************************************
    * @brief H5 Routine to create a main group
    * @param _nameMainGroup Main group name starting with '\'
    * @date 19.03.2020
    * @author Harikrishnan K. Sreekumar
    **********************************************************************/
    void createGroup(std::string _nameMainGroup);
    /**********************************************************************
    * @brief H5 Routine to create a secondary group
    * @param _nameMainGroup Main group name starting with '\'
    * @param _nameSecoGroup Secondary group name starting with '\'
    * @date 19.03.2020
    * @author Harikrishnan K. Sreekumar
    **********************************************************************/
    void createSecoGroup(std::string _nameMainGroup, std::string _nameSecoGroup);
    /**********************************************************************
    * @brief H5 Routine to append an integer vector
    * @param _nameMainGroup Main roup name starting with '\'
    * @param _nameMatGroup Matrix name
    * @param _nameVecGroup Vector name
    *      The structure of resulting H5 :
    *           H5 RootFile:
    *                   -- _nameMainGroup
    *                          -- _nameMatGroup 
    *                                   -- _nameVecGroup
    * @date 19.03.2020
    * @author Harikrishnan K. Sreekumar
    **********************************************************************/
    void appendIntegerVector(std::vector<int>& _vector, std::string _nameMainGroup,  std::string _nameMatGroup, std::string _nameVecGroup);
    /**********************************************************************
    * @brief H5 Routine to append an complex vector
    * @param _vector Vector to be appendend to HDF5 file
    * @param _nameMainGroup Main roup name starting with '\'
    * @param _nameMatGroup Matrix name
    * @param _nameVecGroup Vector name
    *      The structure of resulting H5 :
    *           H5 RootFile:
    *                   -- _nameMainGroup
    *                          -- _nameMatGroup
    *                                   -- _nameVecGroup
    * @date 19.03.2020
    * @author Harikrishnan K. Sreekumar
    **********************************************************************/
#ifdef HAVE_PETSC
    void appendComplexVector(std::vector<PetscComplex>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup);
#endif
    /**********************************************************************
    * @brief H5 Routine to append a complex matrix
    * @param _vector Vector of matrix [Row-major]
    * @param _m Number of rows
    * @param _n Number of columns
    * @param _nameMainGroup Main roup name starting with '\'
    * @param _nameMatGroup Matrix name
    * @param _nameVecGroup Matrix name
    *      The structure of resulting H5 :
    *           H5 RootFile:
    *                   -- _nameMainGroup
    *                          -- _nameMatGroup
    *                                   -- _nameVecGroup
    * @date 25.08.2020
    * @author Harikrishnan K. Sreekumar
    **********************************************************************/
#ifdef HAVE_PETSC
    void appendComplexMatrix(std::vector<PetscComplex>& _vectorMatrix, size_t _m, size_t _n, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup);
#endif
    /**********************************************************************
    * @brief H5 Routine to append a dense real matrx
    * @param _mat MAtrix to be exported
    * @param _m Number of rows
    * @param _n Number of columns
    * @param _nameMainGroup Main roup name starting with '\'
    * @param _nameMatGroup Matrix name
    * @param _nameVecGroup Vector name
    *      The structure of resulting H5 :
    *           H5 RootFile:
    *                   -- _nameMainGroup
    *                          -- _nameMatGroup
    *                                   -- _nameVecGroup
    * @date 25.04.2020
    * @author Harikrishnan K. Sreekumar
    **********************************************************************/
#ifdef HAVE_PETSC
    void appendRealDenseMatrix(std::vector<PetscReal>& _matrix, int _m, int _n, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup);
#endif
    void appendIntegerMatrix(std::vector<int>& _vectorMatrix, size_t _m, size_t _n, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup);
	void appendDoubleVector(std::vector<double>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup);
    /**********************************************************************
    * @brief Constructor
    * @date 19.03.2020
    * @author Harikrishnan K. Sreekumar
    **********************************************************************/
    cOutputH5();
    /**********************************************************************
    * @brief Destructor
    * @date 19.03.2020
    * @author Harikrishnan K. Sreekumar
    **********************************************************************/
    ~cOutputH5();
    /**********************************************************************
    * @brief Function to return the singleton instance
    * @date 24.03.2020
    * @author Harikrishnan K. Sreekumar
    **********************************************************************/
    static cOutputH5* getInstance() {
            if (!mySingleton) {
                mySingleton = new cOutputH5();
            }
            return mySingleton;
    }

    /**********************************************************************
    * @brief H5 Routine to write a string out of an attribute
    * @param _name Name of attribute
    * @param _path Group path to find the attribute
    * @param _value Value of attribute
    * @date 25.03.2020
    * @author Harikrishnan K. Sreekumar
    **********************************************************************/
    void writeStringAttribute(std::string _name, std::string _path, std::string _value);
private:
#ifdef USE_HDF5
    /// H5 Group
    H5::Group* myGroup;
    /// H5 DataSet
    H5::DataSet* myDataSet;
    /// H5 DataSpace
    H5::DataSpace* myDataSpace;
#endif
    /// Singleton
    static cOutputH5* mySingleton;
    /// Complex Struct
    typedef struct elpasoComplexDouble {
        double real;
        double imag;
    };
};
#endif