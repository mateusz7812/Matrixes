#include "CVct.h"

MyAlgebra::CVct::CVct(int size)
{
	this->size = size;
	vector = new FPTYPE[size];
}

FPTYPE& MyAlgebra::CVct::operator[](int index) const
{
	return vector[index];
}

const MyAlgebra::CVct& MyAlgebra::CVct::operator=(const CVct& rhs)
{
	return *this;
}

const MyAlgebra::CVct& MyAlgebra::CVct::operator=(FPTYPE val)
{
	return *this;
}

MyAlgebra::CVct MyAlgebra::CVct::operator+(const CVct& rhs)
{
	return *this;
}

MyAlgebra::CVct MyAlgebra::CVct::operator-(const CVct& rhs)
{
	return *this;
}

MyAlgebra::CVct MyAlgebra::CVct::operator*(const CMtx& rhs)
{
	return *this;
}

MyAlgebra::CVct MyAlgebra::CVct::operator~()
{
	return *this;
}
