#include "outputh5.h"
#include<iostream>

#define RANK_VEC   1
#define RANK_MAT   2

cOutputH5* cOutputH5::mySingleton = NULL;


cOutputH5::cOutputH5()
{
    H5::Exception::dontPrint();
}

cOutputH5::~cOutputH5()
{
#ifdef USE_HDF5
    myGroup->close();
    delete cFileH5::getFileInstance();
    delete myGroup;
    delete myDataSet;
    delete myDataSpace;
#endif
}

void cOutputH5::createGroup(std::string _nameMainGroup)
{
#ifdef USE_HDF5
    try {
        myGroup = new H5::Group(cFileH5::getFileInstance()->openGroup(_nameMainGroup));
    }
    catch (...) {
        myGroup = new H5::Group(cFileH5::getFileInstance()->createGroup(_nameMainGroup));
    }
#endif
}

void cOutputH5::createSecoGroup(std::string _nameMainGroup, std::string _nameSecoGroup)
{
#ifdef USE_HDF5
    try {
        myGroup->openGroup(_nameMainGroup + _nameSecoGroup);
    }
    catch (...) {
        myGroup->createGroup(_nameMainGroup + _nameSecoGroup);
    }
#endif
}

void cOutputH5::appendIntegerVector(std::vector<int>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef USE_HDF5
        hsize_t dimsf[2];
        dimsf[0] = _vector.size();
        myDataSpace = new H5::DataSpace(RANK_VEC, dimsf);

        try {
            myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup));
        }
        catch (...) {
            myDataSet = new H5::DataSet(cFileH5::getFileInstance()->createDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup, H5::PredType::NATIVE_INT, *myDataSpace));
        }
        myDataSet->write(_vector.data(), H5::PredType::NATIVE_INT);

        myDataSpace->close();
        myDataSet->close();
        delete myDataSet;
        delete myDataSpace;
#endif
}

#ifdef HAVE_PETSC
void cOutputH5::appendComplexVector(std::vector<PetscComplex>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef USE_HDF5
    hsize_t dimsf[2];
    dimsf[0] = _vector.size();
    myDataSpace = new H5::DataSpace(RANK_VEC, dimsf);
    H5::CompType myComplexCompound((size_t)sizeof(elpasoComplexDouble));

    std::vector<elpasoComplexDouble> convertPetscComplex;
    convertPetscComplex.resize(_vector.size());
    for (size_t i = 0; i < _vector.size(); i++)
        convertPetscComplex[i] = {_vector[i].real(), _vector[i].imag()};

    myComplexCompound.insertMember("real", HOFFSET(elpasoComplexDouble, real), H5::PredType::NATIVE_DOUBLE);
    myComplexCompound.insertMember("imag", HOFFSET(elpasoComplexDouble, imag), H5::PredType::NATIVE_DOUBLE);

    try {
        myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup));
    }
    catch (...) {
        myDataSet = new H5::DataSet(cFileH5::getFileInstance()->createDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup, myComplexCompound, *myDataSpace));
    }
    myDataSet->write(convertPetscComplex.data(), myComplexCompound);

    myDataSpace->close();
    myDataSet->close();
    delete myDataSet;
    delete myDataSpace;
#endif
}

void cOutputH5::appendComplexMatrix(std::vector<PetscComplex>& _vectorMatrix, size_t _m, size_t _n, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef USE_HDF5
    hsize_t dimsf[2];
    dimsf[0] = _m;
    dimsf[1] = _n;
    myDataSpace = new H5::DataSpace(RANK_VEC, dimsf);
    H5::CompType myComplexCompound((size_t)sizeof(elpasoComplexDouble));

    std::vector<elpasoComplexDouble> convertPetscComplex;
    convertPetscComplex.resize(_vectorMatrix.size());
    for (size_t i = 0; i < _vectorMatrix.size(); i++)
        convertPetscComplex[i] = { _vectorMatrix[i].real(), _vectorMatrix[i].imag() };

    myComplexCompound.insertMember("real", HOFFSET(elpasoComplexDouble, real), H5::PredType::NATIVE_DOUBLE);
    myComplexCompound.insertMember("imag", HOFFSET(elpasoComplexDouble, imag), H5::PredType::NATIVE_DOUBLE);

    try {
        myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup));
    }
    catch (...) {
        myDataSet = new H5::DataSet(cFileH5::getFileInstance()->createDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup, myComplexCompound, *myDataSpace));
    }
    myDataSet->write(convertPetscComplex.data(), myComplexCompound);

    myDataSpace->close();
    myDataSet->close();
    delete myDataSet;
    delete myDataSpace;
#endif
}
#endif

void cOutputH5::writeStringAttribute(std::string _name, std::string _path, std::string _value)
{
#ifdef USE_HDF5
    H5::DataSpace l_att_dataspace = H5::DataSpace(H5S_SCALAR);

    // Create new string datatype for attribute
    H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters
    
    createGroup(_path);
    H5::Attribute l_attr = myGroup->createAttribute(_name, strdatatype, l_att_dataspace);

    l_attr.write(strdatatype, _value);
#endif
}

#ifdef HAVE_PETSC
void cOutputH5::appendRealDenseMatrix(std::vector<PetscReal>& _matrix, int _m, int _n, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef USE_HDF5
    hsize_t dimsf[2];
    dimsf[0] = _m;
    dimsf[1] = _n;
    myDataSpace = new H5::DataSpace(RANK_MAT, dimsf);
    try {
        myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup));
    }
    catch (...) {
        myDataSet = new H5::DataSet(cFileH5::getFileInstance()->createDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup, H5::PredType::NATIVE_DOUBLE, *myDataSpace));
    }
    myDataSet->write(_matrix.data(), H5::PredType::NATIVE_DOUBLE);

    myDataSpace->close();
    myDataSet->close();
    delete myDataSet;
    delete myDataSpace;
#endif
}
#endif

void cOutputH5::appendIntegerMatrix(std::vector<int>& _vectorMatrix, size_t _m, size_t _n, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef USE_HDF5
    hsize_t dimsf[2];
    dimsf[0] = _m;
    dimsf[1] = _n;
    myDataSpace = new H5::DataSpace(RANK_MAT, dimsf);
    try {
        myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup));
    }
    catch (...) {
        myDataSet = new H5::DataSet(cFileH5::getFileInstance()->createDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup, H5::PredType::NATIVE_INT, *myDataSpace));
    }
    myDataSet->write(_vectorMatrix.data(), H5::PredType::NATIVE_INT);

    myDataSpace->close();
    myDataSet->close();
    delete myDataSet;
    delete myDataSpace;
#endif
}

void cOutputH5::appendDoubleVector(std::vector<double>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef USE_HDF5
    hsize_t dimsf[2];
    dimsf[0] = _vector.size();
    myDataSpace = new H5::DataSpace(RANK_VEC, dimsf);

    try {
        myDataSet = new H5::DataSet(cFileH5::getFileInstance()->openDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup));
    }
    catch (...) {
        myDataSet = new H5::DataSet(cFileH5::getFileInstance()->createDataSet(_nameMainGroup + _nameMatGroup + _nameVecGroup, H5::PredType::NATIVE_DOUBLE, *myDataSpace));
    }
    myDataSet->write(_vector.data(), H5::PredType::NATIVE_DOUBLE);

    myDataSpace->close();
    myDataSet->close();
    delete myDataSet;
    delete myDataSpace;
#endif
}