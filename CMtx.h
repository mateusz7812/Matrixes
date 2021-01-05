#pragma once
#include <string>

namespace UnitTestMatrixes {
	class UnitTestMatrixes;
}

namespace MyAlgebra2
{

	class CMtx
	{
	private:
		CMtx& matrix = *this;
		float** row_ptr;
		int     row_cnt;
		int     col_cnt;
		
	
	public:
		static const float ALG_PRECISION;

		// =========================================================================
		// KONSTRUKTORY:
		// =========================================================================
		
		// Tworzy macierz z mo�liwo�ci� losowej inicjalizacji
		CMtx(int row_cnt, int col_cnt, bool rand_init = false);

		// Tworzy kwadratow� macierz diagonaln�
		CMtx(int row_cnt, float diagonal);

		CMtx(const CMtx& rhs);

		CMtx(CMtx&& rhs);

		// Je�li potrzeba - nale�y zadeklarowa� i zaimplementowa� inne konstruktory
		~CMtx();


		// =========================================================================
		// OPERATORY PRZYPISANIA:
		// =========================================================================

		const CMtx& operator=(const CMtx& rhs);

		// Zamiana macierzy na macierz diagonaln� 
		const CMtx& operator=(float diagonal);

		// Operator pzzenosz�cy
		const CMtx& operator=(CMtx&& rhs);


		// =========================================================================
		// INDEKSOWANIE MACIERZY
		// =========================================================================

		float*& operator[](int row_ind) const; //TODO: remove & , const

		// =========================================================================
		// OPERACJE ALGEBRAICZNE
		// =========================================================================

		CMtx operator*(const CMtx& rhs) const;

		// Mno�enie macierzy przez sta��
		CMtx operator*(float multiplier) const;

		CMtx operator+(const CMtx& rhs) const;
		CMtx operator-(const CMtx& rhs) const;

		// Minus unarny - zmiana znaku wszystkich wsp�czynnik�w macierzy
		CMtx operator-() const;

		// Transponowanie macierzy
		CMtx operator~() const;

		// Akceptuje tylko power >= 0:
		//    power = 0  - zwraca macierz jednostkow�
		//    power = 1  - zwraca kopi� macierzy
		//    power > 1  - zwraca iloczyn macierzy 
		CMtx operator^(int power) const;

		// Por�wnywanie macierzy z dok�adno�ci� do sta�ej ALG_PRECISION
		bool operator==(const CMtx& rhs) const;

		// Tylko do cel�w testowych - wypisuje macierz wierszami na stdout
		void display() const;

		friend CMtx operator*(float multiplier, const MyAlgebra2::CMtx &rhs );
		friend float det(const CMtx& rhs);

		CMtx reversed() const;
		CMtx map(const CMtx& rhs, float fun(float, float)) const;
		void move(CMtx& rhs);
		void copy(const CMtx& rhs);
		float** create_matrix(int row_cnt, int col_cnt);
		float** create_matrix(int row_cnt, int col_cnt, bool rand_init);
		void make_diagonal(int row_cnt, int col_cnt, float diagonal, float** ptr);
		float** create_matrix(int row_cnt, int col_cnt, float diagonal);
		std::string to_string();
		void change_rows(int first, int second);
		friend class UnitTestMatrixes::UnitTestMatrixes;
	};

	CMtx operator*(float multiplier, const CMtx& rhs);
	float det(const CMtx& rhs);
	CMtx** create_sub_matrixes_for_determinant(const CMtx& rhs, const int size);
}

