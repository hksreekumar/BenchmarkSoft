#include "fileh5.h"
#include<iostream>

cFileH5* cFileH5::mySingleton = NULL;
#ifdef USE_HDF5
H5::H5File* cFileH5::myH5File = NULL;
#endif // USE_HDF5
std::string cFileH5::myFileName = "";

cFileH5::cFileH5()
{
	// Empty
}

cFileH5::~cFileH5()
{
}

void cFileH5::openContainer(ELPASO_H5MODE _mode)
{
	switch (_mode)
	{
	case ELPASO_H5_READONLY:
#ifdef USE_HDF5
		// Existing file is opened with read-only access. If file
		// does not exist, H5Fopen fails.
		myH5File = new H5::H5File(myFileName, H5F_ACC_RDONLY);
#endif
		break;
	case ELPASO_H5_READWRITE:
#ifdef USE_HDF5
		// Existing file is opened with read - write access. If file
		// does not exist, H5Fopen fails.
		myH5File = new H5::H5File(myFileName, H5F_ACC_RDWR);
#endif
		break;
	case ELPASO_H5_READWRITE_FORCE:
#ifdef USE_HDF5
		// If file already exists, file is opened with read-write
		// access and new data overwrites existing data,
		// destroying all prior content, i.e., file content is truncated
		// upon opening. If file does not exist, it is created and
		// opened with read-write access.
		myH5File = new H5::H5File(myFileName, H5F_ACC_TRUNC);
#endif
		break;
	default:
		std::cerr << "	! Invalid H5 open mode !" << std::endl;
		break;
	}
}

void cFileH5::closeContainer()
{
#ifdef USE_HDF5
	myH5File->close();
#endif
}




