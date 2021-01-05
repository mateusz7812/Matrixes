#pragma once

typedef float FPTYPE;

#include "CMtx.h"

namespace MyAlgebra
{
	class CMtx;

	class CVct
	{
	private:
		FPTYPE* vector;
		int     size;

	public:
		CVct(int size);

		FPTYPE& operator[](int index) const;

		const CVct& operator=(const CVct& rhs);
		const CVct& operator=(FPTYPE val);

		CVct operator+(const CVct& rhs);
		CVct operator-(const CVct& rhs);

		CVct operator*(const CMtx& rhs);

		// Transpozycja - zamiana wektora wierszowego na kolumnowy i odwrotnie
		CVct operator~();
	};
	
}
