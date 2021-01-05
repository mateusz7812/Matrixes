#pragma once
#include "../MyAlgebra/CMtx.h"
#include <string>

namespace UnitTestMatrixes {
	class UnitTestMatrixes;
}

typedef float FPTYPE;

namespace RefAlgebra
{
	class CVctRef;

	class CMtxRef
	{
	private:
		CMtxRef& matrix = *this;
		FPTYPE** row_ptr;
		int     row_cnt;
		int     col_cnt;

		friend CVctRef;
	
	public:
		static const FPTYPE ALG_PRECISION ;

		// =========================================================================
		// KONSTRUKTORY:
		// =========================================================================

		// Tworzy pustą macierz 0x0
		CMtxRef();
		
		// Tworzy macierz z mo�liwo�ci� losowej inicjalizacji
		CMtxRef(int row_cnt, int col_cnt, bool rand_init = false);

		// Tworzy kwadratow� macierz diagonaln�
		CMtxRef(int row_cnt, FPTYPE diagonal);

		CMtxRef(const CMtxRef& rhs);

		CMtxRef(CMtxRef&& rhs);

		// Je�li potrzeba - nale�y zadeklarowa� i zaimplementowa� inne konstruktory
		~CMtxRef();


		// =========================================================================
		// OPERATORY PRZYPISANIA:
		// =========================================================================

		const CMtxRef& operator=(const CMtxRef& rhs);

		// Zamiana macierzy na macierz diagonaln� 
		const CMtxRef& operator=(const FPTYPE diagonal);

		// Operator pzzenosz�cy
		const CMtxRef& operator=(CMtxRef&& rhs);


		// =========================================================================
		// INDEKSOWANIE MACIERZY
		// =========================================================================

		FPTYPE*& operator[](int row_ind) const;

		// =========================================================================
		// OPERACJE ALGEBRAICZNE
		// =========================================================================

		// Mnozenie macierzy przez wektor, rhs musi by� wektorem kolumnowym
		CVctRef operator*(const CVctRef& rhs) const;

		CMtxRef operator*(const CMtxRef& rhs) const;

		// Mno�enie macierzy przez sta��
		CMtxRef operator*(FPTYPE multiplier) const;

		CMtxRef operator+(const CMtxRef& rhs) const;
		CMtxRef operator-(const CMtxRef& rhs) const;

		// Minus unarny - zmiana znaku wszystkich wsp�czynnik�w macierzy
		CMtxRef operator-() const;

		// Transponowanie macierzy
		CMtxRef operator~() const;

		// Akceptuje tylko power >= -1:
		//    power = -1 - zwraca macierz odwr�con�
		//    power = 0  - zwraca macierz jednostkow�
		//    power = 1  - zwraca kopi� macierzy
		//    power > 1  - zwraca iloczyn macierzy 
		CMtxRef operator^(int power) const;

		CMtxRef reversed() const;

		// Por�wnywanie macierzy z dok�adno�ci� do sta�ej ALG_PRECISION
		bool operator==(const CMtxRef& rhs) const;

		// Tylko do cel�w testowych - wypisuje macierz wierszami na stdout
		void display() const;

		friend CMtxRef operator*( FPTYPE multiplier, const RefAlgebra::CMtxRef &rhs );
		friend float det(const CMtxRef& rhs);

		CMtxRef map(const CMtxRef& rhs, FPTYPE fun(FPTYPE, FPTYPE)) const;
		void move(CMtxRef& rhs);
		void copy(const CMtxRef& rhs);
		FPTYPE** create_matrix(int row_cnt, int col_cnt);
		FPTYPE** create_matrix(int row_cnt, int col_cnt, bool rand_init);
		void make_diagonal(int row_cnt, int col_cnt, FPTYPE diagonal, FPTYPE** ptr);
		FPTYPE** create_matrix(int row_cnt, int col_cnt, FPTYPE diagonal);
		std::string to_string();
		void change_rows(int first, int second);
		friend class UnitTestMatrixes::UnitTestMatrixes;
	};

	CMtxRef operator*(FPTYPE multiplier, const CMtxRef& rhs);
	float det(const CMtxRef& rhs);
	CMtxRef** create_sub_matrixes_for_determinant(const CMtxRef& rhs, const int size);
}

