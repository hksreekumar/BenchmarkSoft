#ifndef INFAM_INPUTH5_H
#define INFAM_INPUTH5_H

#include <vector>
#include "fileh5.h"

#ifdef HAVE_PETSC
#include "petscconf.h"
#include "petscksp.h"
#include "petscao.h"
#include "slepceps.h"
#endif // HAVE_PETSC


#ifdef USE_INTEL_MKL
#include "mkl_pardiso.h"
#include "mkl_types.h"
#endif

/**
* @brief Singleton H5 class supporting H5 read
* @author Harikrishnan Sreekumar
* @date 23.03.2020
*
* cOutputH5 reads data from a h5 file.
*/

class cInputH5: public cFileH5
{
public:
	/**********************************************************************
	* @brief Constructor
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	cInputH5();
	/**********************************************************************
	* @brief Destructor
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	~cInputH5();
	/**********************************************************************
	* @brief H5 Routine to read an integer vector
	* @param _vector Vector handle
	* @param _nameMainGroup Main roup name starting with '\'
	* @param _nameMatGroup Matrix name
	* @param _nameVecGroup Vector name
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	void readIntegerVector(std::vector<int>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup);
	/**********************************************************************
	* @brief H5 Routine to read a dense integer matrix
	* @param _vector Vector handle containing all the matrix details
	* @param _m Number of rows read
	* @param _n Number of columns read
	* @param _nameMainGroup Main roup name starting with '\'
	* @param _nameMatGroup Matrix name
	* @param _nameVecGroup Vector name
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	void readDenseDoubleMatrix(std::vector<double>& _vector, int& _m, int& _n, std::string _nameFirstLevel, std::string _nameSecondLevel, std::string _nameThirdLevel);
	/**********************************************************************
	* @brief H5 Routine to read an double vector
	* @param _vector Vector handle
	* @param _nameMainGroup Main roup name starting with '\'
	* @param _nameMatGroup Matrix name
	* @param _nameVecGroup Vector name
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	void readDoubleVector(std::vector<double>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup);
	/**********************************************************************
	* @brief H5 Routine to read an complex vector
	* @param _vector Vector handle
	* @param _nameMainGroup Main group name starting with '\'
	* @param _nameMatGroup Matrix name
	* @param _nameVecGroup Vector name
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
#ifdef HAVE_PETSC
	void readComplexVector(std::vector<PetscComplex>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup);
#endif
	/**********************************************************************
	* @brief H5 Routine to read an complex vector using intel mkl
	* @param _vector Vector handle
	* @param _nameMainGroup Main group name starting with '\'
	* @param _nameMatGroup Matrix name
	* @param _nameVecGroup Vector name
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
#ifdef USE_INTEL_MKL
	void readComplexVector(std::vector<MKL_Complex16>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup);
#endif
	/**********************************************************************
	* @brief Function to return the singleton instance
	* @date 19.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	static cInputH5* getInstance() {
		if (!mySingleton) {
			mySingleton = new cInputH5();
		}
		return mySingleton;
	}
	/**********************************************************************
	* @brief H5 Routine to read a integer out of an attribute contained in a group
	* @param _name Name of attribute
	* @param _path Group path to find the attribute
	* @return Read int
	* @date 25.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	int readIntegerAttributeFromGroup(std::string _name, std::string _path);
	/**********************************************************************
	* @brief H5 Routine to read a string out of an attribute contained in a group
	* @param _name Name of attribute
	* @param _path Group path to find the attribute
	* @return Read string
	* @date 25.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	std::string readStringAttributeFromGroup(std::string _name, std::string _path);
	/**********************************************************************
	* @brief H5 Routine to read a string out of an attribute contained in a dataset
	* @param _name Name of attribute
	* @param _path Group path to find the attribute
	* @param _dataset Name of dataset
	* @return Read string
	* @date 03.07.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	std::string readStringAttributeFromDataset(std::string _name, std::string _path, std::string _dataset);

	/**********************************************************************
	* @brief H5 Routine to read a integer out of an attribute contained in a dataset
	* @param _name Name of attribute
	* @param _path Group path to find the attribute
	* @param _dataset Name of dataset
	* @return Read integer
	* @date 20.10.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	int readIntegerAttributeFromDataset(std::string _name, std::string _path, std::string _dataset);

	/**********************************************************************
	* @brief H5 Routine to get number of members of type datasets in given group path
	* @param _path Group path
	* @return number of members
	* @date 20.10.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	int getNumberOfMembersInGroup(std::string _path);
private:
	/// Singleton
	static cInputH5* mySingleton;
	/// H5 DataSet
	H5::DataSet* myDataSet;
	/// H5 DataSpace
	H5::DataSpace* myDataSpace;
	/// Complex Struct
	typedef struct elpasoComplexDouble {
		double real;
		double imag;
	};
};
#endif // !INFAM_INPUTH5_H

