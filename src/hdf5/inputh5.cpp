#include "inputh5.h"
#include <string>
#include <iostream>

#define RANK_VEC   1

cInputH5* cInputH5::mySingleton = NULL;

cInputH5::cInputH5()
{
}

cInputH5::~cInputH5()
{
}

void cInputH5::readIntegerVector(std::vector<int>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef USE_HDF5
    myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup));

    myDataSpace = new H5::DataSpace(myDataSet->getSpace());
    int rank = myDataSpace->getSimpleExtentNdims();
    hsize_t dims_out[2];
    int ndims = myDataSpace->getSimpleExtentDims(dims_out, NULL);
    _vector.resize(dims_out[0]);
    myDataSet->read(&_vector[0], H5::PredType::NATIVE_INT);

    myDataSpace->close();
    myDataSet->close();
    delete myDataSet;
    delete myDataSpace;
#endif
}

void cInputH5::readDenseDoubleMatrix(std::vector<double>& _vector, int& _m, int& _n, std::string _nameFirstLevel, std::string _nameSecondLevel, std::string _nameThirdLevel)
{
#ifdef USE_HDF5
    myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_nameFirstLevel + _nameSecondLevel + _nameThirdLevel));

    myDataSpace = new H5::DataSpace(myDataSet->getSpace());
    int rank = myDataSpace->getSimpleExtentNdims();
    hsize_t dims_out[2];
    int ndims = myDataSpace->getSimpleExtentDims(dims_out, NULL);
    _vector.resize(dims_out[0]* dims_out[1]);

    myDataSet->read(&_vector[0], H5::PredType::NATIVE_DOUBLE);

    _m = (int)dims_out[0];
    _n = (int)dims_out[1];

    myDataSpace->close();
    myDataSet->close();
    delete myDataSet;
    delete myDataSpace;
#endif 
}

void cInputH5::readDoubleVector(std::vector<double>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef USE_HDF5
    myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup));

    myDataSpace = new H5::DataSpace(myDataSet->getSpace());
    int rank = myDataSpace->getSimpleExtentNdims();
    hsize_t dims_out[2];
    int ndims = myDataSpace->getSimpleExtentDims(dims_out, NULL);
    _vector.resize(dims_out[0]);
    myDataSet->read(&_vector[0], H5::PredType::NATIVE_DOUBLE);

    myDataSpace->close();
    myDataSet->close();
    delete myDataSet;
    delete myDataSpace;
#endif
}

#ifdef HAVE_PETSC
void cInputH5::readComplexVector(std::vector<PetscComplex>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef USE_HDF5
    myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup));

    myDataSpace = new H5::DataSpace(myDataSet->getSpace());
    int rank = myDataSpace->getSimpleExtentNdims();
    hsize_t dims_out[2];
    int ndims = myDataSpace->getSimpleExtentDims(dims_out, NULL);

    H5::CompType myComplexCompound((size_t)sizeof(elpasoComplexDouble));

    std::vector<elpasoComplexDouble> convertPetscComplex;
    convertPetscComplex.resize(dims_out[0]);
    _vector.resize(dims_out[0]);

    myComplexCompound.insertMember("real", HOFFSET(elpasoComplexDouble, real), H5::PredType::NATIVE_DOUBLE);
    myComplexCompound.insertMember("imag", HOFFSET(elpasoComplexDouble, imag), H5::PredType::NATIVE_DOUBLE);
    
    myDataSet->read(&convertPetscComplex[0], myComplexCompound);
    for (size_t i = 0; i < convertPetscComplex.size(); i++)
        _vector[i] = { convertPetscComplex[i].real, convertPetscComplex[i].imag };

    myDataSpace->close();
    myDataSet->close();
    delete myDataSet;
    delete myDataSpace;
#endif
}

int cInputH5::readIntegerAttributeFromGroup(std::string _name, std::string _path)
{
#ifdef USE_HDF5
try {
    int readBuffer = 0;
    H5::Group myGroup = H5::Group(cFileH5::getFileInstance()->openGroup(_path));

    H5::Attribute l_attr = myGroup.openAttribute(_name);
    H5::DataType type = l_attr.getDataType();
    l_attr.read(type, &readBuffer);

    return readBuffer;
}
catch (H5::AttributeIException attribute_not_found) {
    std::cout << "Attribute not found" << std::endl;
    return -1000;
}
#endif
}

std::string cInputH5::readStringAttributeFromGroup(std::string _name, std::string _path)
{
#ifdef USE_HDF5
try {
    H5::Group myGroup = H5::Group(cFileH5::getFileInstance()->openGroup(_path));

    H5::Attribute l_attr = myGroup.openAttribute(_name);
    H5::DataType type = l_attr.getDataType();
    H5std_string strreadbuf("");
    l_attr.read(type, strreadbuf);

    return strreadbuf;
}
catch (H5::AttributeIException attribute_not_found) {
    std::cout << "Attribute not found" << std::endl;
    return "Not found";
}
#endif
}

std::string cInputH5::readStringAttributeFromDataset(std::string _name, std::string _path, std::string _dataset)
{
#ifdef USE_HDF5
    myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_path + _dataset));
    H5::Attribute myatt_out = myDataSet->openAttribute(_name);
    H5::DataType type = myatt_out.getDataType();
    H5std_string strreadbuf("");
    myatt_out.read(type, strreadbuf);

    delete myDataSet;

    return strreadbuf;
#endif
}

int cInputH5::readIntegerAttributeFromDataset(std::string _name, std::string _path, std::string _dataset)
{
#ifdef USE_HDF5
    int readBuffer = 0;
    myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_path + _dataset));

    H5::Attribute l_attr = myDataSet->openAttribute(_name);
    H5::DataType type = l_attr.getDataType();
    l_attr.read(type, &readBuffer);

    delete myDataSet;

    return readBuffer;
#endif
}

int cInputH5::getNumberOfMembersInGroup(std::string _path)
{
#ifdef USE_HDF5
    try {
        H5::Group myGroup = H5::Group(cFileH5::getFileInstance()->openGroup(_path));
        return myGroup.getNumObjs();
    }
    catch (H5::GroupIException group_not_found) {
        std::cout << "Group not found" << std::endl;
        return -1;
    }
#endif
}
#endif

#ifdef USE_INTEL_MKL
void cInputH5::readComplexVector(std::vector<MKL_Complex16>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef USE_HDF5
    myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup));

    myDataSpace = new H5::DataSpace(myDataSet->getSpace());
    int rank = myDataSpace->getSimpleExtentNdims();
    hsize_t dims_out[2];
    int ndims = myDataSpace->getSimpleExtentDims(dims_out, NULL);

    H5::CompType myComplexCompound((size_t)sizeof(elpasoComplexDouble));

    std::vector<elpasoComplexDouble> convertPetscComplex;
    convertPetscComplex.resize(dims_out[0]);
    _vector.resize(dims_out[0]);

    myComplexCompound.insertMember("real", HOFFSET(elpasoComplexDouble, real), H5::PredType::NATIVE_DOUBLE);
    myComplexCompound.insertMember("imag", HOFFSET(elpasoComplexDouble, imag), H5::PredType::NATIVE_DOUBLE);

    myDataSet->read(&convertPetscComplex[0], myComplexCompound);
    for (size_t i = 0; i < convertPetscComplex.size(); i++)
        _vector[i] = { convertPetscComplex[i].real, convertPetscComplex[i].imag };

    myDataSpace->close();
    myDataSet->close();
    delete myDataSet;
    delete myDataSpace;
#endif
}
#ifdef