
#include <cassert>
#include <iostream>


#include "MyAlgebra/CMtx.h"
#include "RefAlgebra/CMtxRef.h"

#include <stdlib.h>
#include <stdio.h>


// ===================================================================
// FUNKCJE DO POMIARU CZASU
// ===================================================================

#ifdef _WIN32
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif
#include <time.h>
#include<math.h>

#include "RefAlgebra/CVctRef.h"


double mygettime(void) {
# ifdef _WIN32
	struct _timeb tb;
	_ftime64_s(&tb);
	return static_cast<double>(tb.time) + 0.001 * static_cast<double>(tb.millitm);
# else
	struct timeval tv;
	if (gettimeofday(&tv, 0) < 0) {
		perror("oops");
	}
	return (double)tv.tv_sec + (0.000001 * (double)tv.tv_usec);
# endif
}

// ===================================================================
// FUNKCJA OCENY CZASU WYKONANIA
// ===================================================================

// Definiujemy szablon aby �atwiej uruchamia� testy dla r�znych implementacji
// klasy. R�ne implementacje b�d� umieszczone w r�nych przestrzeniach nazw.
template<typename T>
double test()
{
	// Przyk�adowe testowe obliczenie macierzowe. Podobne obliczenia b�d� 
	// u�ywane do oceny efektywno�ci implementacji w konkursie.
	const int SIZE = 100;
	const int ITER_CNT = 10;

	T A(SIZE, SIZE, true);
	T B(SIZE, SIZE, true);
	T W(1, 1);
	double t1 = mygettime();

	for (int i = 0; i < ITER_CNT; i++)
	{
		//B = ((/*0.1 */ i) * A + B * B); //* 1.e-4;
		B = -B * ~(A + B);
	}
	W = (B - A);

	double exec_time = mygettime() - t1;

	W.display();

	return exec_time;
}

void _tmain()
{
	double t_prog = test<MyAlgebra::CMtx>();

	double t_ref = test<RefAlgebra::CMtxRef>();

	printf("Czas wykonania referencyjny: %7.2lfs\n", t_ref);
	printf("Czas wykonania testowany:    %7.2lfs\n", t_prog);
	printf("Wsp�czynnik przyspieszenia Q: %5.2lf", t_ref / t_prog);
}


int main()
{
	std::cout << RefAlgebra::CMtxRef::ALG_PRECISION;
    //_tmain();
}
