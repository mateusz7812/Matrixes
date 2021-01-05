#pragma once

typedef float FPTYPE;

#include "CMtxRef.h"

namespace RefAlgebra
{
	class CMtxRef;

	class CVctRef
	{
	public:
		CVctRef(int size);
		CVctRef(const CVctRef& rhs);
		CVctRef(CVctRef&& rhs);

		FPTYPE& operator[](int ind) const;

		const CVctRef& operator=(const CVctRef& rhs);
		const CVctRef& operator=(CVctRef&& rhs);
		const CVctRef& operator=(FPTYPE val);

		CVctRef operator+(const CVctRef& rhs);
		CVctRef operator-(const CVctRef& rhs);

		CVctRef operator*(const CMtxRef& rhs);

		// Transpozycja - zamiana wektora wierszowego na kolumnowy i odwrotnie
		CVctRef operator~();

		bool is_column();
		void copy(const CVctRef& rhs);
		void move(CVctRef& rhs);


	private:
		FPTYPE* vector;
		int     size;
		bool column;
		CVctRef map(const CVctRef& rhs, FPTYPE fun(FPTYPE, FPTYPE));
	};

	
}
