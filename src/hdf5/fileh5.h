
#ifndef INFAM_FILEH5_H
#define INFAM_FILEH5_H

#include<string>

#ifdef USE_HDF5
#include "H5Cpp.h"
#endif // USE_HDF5

enum ELPASO_H5MODE
{
	ELPASO_H5_READONLY,
	ELPASO_H5_READWRITE,
	ELPASO_H5_READWRITE_FORCE
};

/**
* @brief Implements singleton pattern for the global elpaso H5 instance
* @author Harikrishnan Sreekumar
* @date 23.03.2020
*
* cFileH5 is the base class for input and output routines of H5
*/

class cFileH5
{
public:
	/**********************************************************************
	* @brief Constructor
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	cFileH5();
	/**********************************************************************
	* @brief Destructor
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	~cFileH5();
	/**********************************************************************
	* @brief Base H5 Routine to set h5 file name
	* @param _filename Name of file
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	void setGlobalFileName(std::string _filename) { myFileName = _filename; }
	/**********************************************************************
	* @brief Base H5 Routine to get h5 file name
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	std::string getGlobalFileName() { return myFileName; }
	/**********************************************************************
	* @brief Base H5 Routine to close a h5 file
	* @param _mode Mode for opening container [READ-ONLY or READ/WRITE]
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	void openContainer(ELPASO_H5MODE _mode);
	/**********************************************************************
	* @brief Base H5 Routine to close a h5 file
	* @date 23.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	virtual void closeContainer();

	/**********************************************************************
	* @brief Function to return the static filename instance
	* @date 24.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
#ifdef USE_HDF5
	static H5::H5File* getFileInstance() {
		if (myH5File) {
			return myH5File;
		}
		return NULL;
	};
#endif // USE_HDF5
	/**********************************************************************
	* @brief Function to return the singleton instance
	* @date 24.03.2020
	* @author Harikrishnan K. Sreekumar
	**********************************************************************/
	static cFileH5* getInstance() {
		if (!mySingleton) {
			mySingleton = new cFileH5();
		}
		return mySingleton;
	}
private:
	/// Singleton
	static cFileH5* mySingleton;
#ifdef USE_HDF5
	/// H5 Root file
	static H5::H5File* myH5File;
#endif
	/// File name
	static std::string myFileName;
};
#endif