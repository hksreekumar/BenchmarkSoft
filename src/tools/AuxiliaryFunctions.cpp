#include <AuxiliaryFunctions.h>
#include <fstream>
#include <iostream>
#include <limits>

AuxiliaryFunctions::AuxiliaryFunctions()
{
}

AuxiliaryFunctions::~AuxiliaryFunctions()
{
}


std::vector<double> AuxiliaryFunctions::readDoubleVectorDatFormat(std::string _fileName){
    std::cout<<" > AuxFun: Reading file " << _fileName <<"...";
    std::vector<double> temp;
    
    std::ifstream fin (_fileName);
    double num;
    fin >> std::scientific;
    
    while(fin >> num)
        temp.push_back(num);
    std::cout<<" Finished." << std::endl;
    return temp;
    
}

std::vector<int> AuxiliaryFunctions::readIntegerVectorDatFormat(std::string _fileName){
    std::cout<<" > AuxFun: Reading file " << _fileName <<"...";
    std::vector<int> temp;
    
    std::ifstream fin (_fileName);
    double num;
    fin >> std::scientific;

    while(fin >> num)
        temp.push_back((int)num);
    std::cout<<" Finished." << std::endl;
    return temp;
}

void AuxiliaryFunctions::writeDoubleVectorDatFormat(std::string _fileName, std::vector<double> &_vector) {
	std::cout << ">> Writing " << _fileName << "#" << _vector.size() <<"..." << std::endl;
	size_t ii_couter;
	std::ofstream myfile;
	myfile.open(_fileName);
	myfile.precision(std::numeric_limits<double>::digits10 + 1);
	myfile << std::scientific;
	for (ii_couter = 0; ii_couter < _vector.size(); ii_couter++)
	{
		myfile << _vector[ii_couter] << std::endl;
	}
	myfile << std::endl;
	myfile.close();

}

#ifdef USE_INTEL_MKL
std::vector<MKL_Complex16> AuxiliaryFunctions::readComplexVectorDatFormat(std::string _fileName) {
    std::cout<<" > AuxFun: Reading complex file " << _fileName <<"...";
    std::vector<MKL_Complex16> temp;
    
    std::ifstream fin (_fileName);
    double num_r;
    double num_i;
    fin >> std::scientific;

    while(fin >> num_r >> num_i){
        MKL_Complex16 __t; __t.real = num_r; __t.imag = num_i;
        temp.push_back(__t);
    }
    std::cout<<" Finished." << std::endl;
    return temp;
}
#endif

#ifdef USE_PETSC
std::vector<PetscComplex> AuxiliaryFunctions::readPetscDoubleVectorDatFormatAsComplex(std::string _fileName) {
    std::cout<<" > AuxFun: Reading complex file " << _fileName <<"...";
    std::vector<PetscComplex> temp;
    
    std::ifstream fin (_fileName);
    double num_r;
    fin >> std::scientific;

    while(fin >> num_r){
        PetscComplex __t = num_r;
        temp.push_back(__t);
    }
    std::cout<<" Finished." << std::endl;
    return temp;
}

std::vector<PetscComplex> AuxiliaryFunctions::readPetscComplexVectorDatFormat(std::string _fileName) {
    std::cout<<" > AuxFun: Reading complex file " << _fileName <<"...";
    std::vector<PetscComplex> temp;
    
    std::ifstream fin (_fileName);
    double num_r;
    double num_i;
    fin >> std::scientific;

    while(fin >> num_r >> num_i){
        PetscComplex __t = num_r + num_i * PETSC_i;
        temp.push_back(__t);
    }
    std::cout<<" Finished." << std::endl;
    return temp;
}
#endif

#ifdef USE_INTEL_MKL
void AuxiliaryFunctions::writeComplexVectorDatFormat(std::string _fileName, std::vector<MKL_Complex16> &_vector) {
	std::cout << ">> Writing " << _fileName << "#" << _vector.size() <<"..." << std::endl;
	size_t ii_couter;
	std::ofstream myfile;
	myfile.open(_fileName);
	myfile.precision(std::numeric_limits<double>::digits10 + 1);
	myfile << std::scientific;
	for (ii_couter = 0; ii_couter < _vector.size(); ii_couter++)
	{
		myfile << _vector[ii_couter].real << "\t" << _vector[ii_couter].imag << std::endl;
	}
	myfile << std::endl;
	myfile.close();

}
#endif

#ifdef USE_PETSC
void AuxiliaryFunctions::writeComplexVectorDatFormat(std::string _fileName, std::vector<PetscComplex> &_vector) {
	std::cout << ">> Writing " << _fileName << "#" << _vector.size() <<"..." << std::endl;
	size_t ii_couter;
	std::ofstream myfile;
	myfile.open(_fileName);
	myfile.precision(std::numeric_limits<double>::digits10 + 1);
	myfile << std::scientific;
	for (ii_couter = 0; ii_couter < _vector.size(); ii_couter++)
	{
		myfile << _vector[ii_couter].real() << "\t" << _vector[ii_couter].imag() << std::endl;
	}
	myfile << std::endl;
	myfile.close();

}
#endif

#ifdef USE_INTEL_MKL
std::vector<double> AuxiliaryFunctions::convertComplexVectorToDouble(std::vector<MKL_Complex16> &_vector) {
    std::cout<<">> Converting Complex To Double Vector ... ";
    std::vector<double> temp(_vector.size());
    for(int i = 0 ; i< _vector.size(); i++) {
        temp[i]=_vector[i].real;
    }
    std::cout<<" Finished." << std::endl;
    
    return temp;
}
#endif

#ifdef USE_INTEL_MKL
std::vector<MKL_Complex16> AuxiliaryFunctions::convertDoubleVectorToComplex(std::vector<double> &_vector) {
    std::cout<<">> Converting Double To Complex Vector ... ";
    std::vector<MKL_Complex16> temp(_vector.size());
    for(int i = 0 ; i< _vector.size(); i++) {
        temp[i].real =_vector[i];
    }
    std::cout<<" Finished." << std::endl;
    return temp;
}
#endif
