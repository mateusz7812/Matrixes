#pragma once
typedef float FPTYPE;

namespace MyAlgebra
{
	class CVct;

	class CMtx
	{
	private:
		FPTYPE** row_ptr;
		int     row_cnt;
		int     col_cnt;
		
		friend class CVct;

	public:
		static const FPTYPE ALG_PRECISION;

		// =========================================================================
		// KONSTRUKTORY:
		// =========================================================================

		// Tworzy macierz z mo�liwo�ci� losowej inicjalizacji
		CMtx(int row_cnt, int col_cnt, bool rand_init = false);

		// Tworzy kwadratow� macierz diagonaln�
		CMtx(int row_cnt, FPTYPE diagonal);

		CMtx(const CMtx& rhs);

		// Je�li potrzeba - nale�y zadeklarowa� i zaimplementowa� inne konstruktory
		~CMtx();


		// =========================================================================
		// OPERATORY PRZYPISANIA:
		// =========================================================================

		const CMtx& operator=(const CMtx& rhs);

		// Zamiana macierzy na macierz diagonaln� 
		const CMtx& operator=(const FPTYPE diagonal);

		// Operator pzzenosz�cy
		const CMtx& operator=(CMtx&& rhs);


		// =========================================================================
		// INDEKSOWANIE MACIERZY
		// =========================================================================

		FPTYPE* operator[](int row_ind);

		// =========================================================================
		// OPERACJE ALGEBRAICZNE
		// =========================================================================

		// Mnozenie macierzy przez wektor, rhs musi by� wektorem kolumnowym
		CVct operator*(const CVct& rhs) const;

		CMtx operator*(const CMtx& rhs) const;

		// Mno�enie macierzy przez sta��
		CMtx operator*(FPTYPE multiplier) const;

		CMtx operator+(const CMtx& rhs) const;
		CMtx operator-(const CMtx& rhs) const;

		// Minus unarny - zmiana znaku wszystkich wsp�czynnik�w macierzy
		CMtx operator-() const;

		// Transponowanie macierzy
		CMtx operator~() const;

		// Akceptuje tylko power >= -1:
		//    power = -1 - zwraca macierz odwr�con�
		//    power = 0  - zwraca macierz jednostkow�
		//    power = 1  - zwraca kopi� macierzy
		//    power > 1  - zwraca iloczyn macierzy 
		CMtx operator^(int power) const;

		// Por�wnywanie macierzy z dok�adno�ci� do sta�ej ALG_PRECISION
		bool operator==(const CMtx&& rhs) const;

		// Tylko do cel�w testowych - wypisuje macierz wierszami na stdout
		void display() const;

		// friend CMtx operator*( FPTYPE multiplier, const CMtx &rhs );
	};

	CMtx operator*(FPTYPE multiplier, const CMtx& rhs);
}

