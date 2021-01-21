#include <iostream>
#include <map>
#include <conio.h>
#include "../matrix_default/CMtx.h"
using namespace MyOptional;
using namespace MyAlgebra;
using namespace std;

typedef int TYPE;

class CProgram
{
	map<string, CMtx<TYPE>*> _matrixes;

	string get_message_for_code(int code)
	{
		switch (code)
		{
		case 0:
			return "Brak błędów";
		case 1:
			return "podana zostałą macierz o niewłaściwych wymiarach";
		case 2:
			return "Podana została zła wartość";
		case 3:
			return "Podana macierz jest osobliwa";
		case 4:
			return "Pliku nie znaleziono";
		default:
			return "Nierozpoznany kod";
		}
	}

	string to_string_all_matrixes()
	{
		stringstream stringstream;
		for (map<string, CMtx<TYPE>*>::value_type pair : _matrixes)
		{
			stringstream << pair.first << ":\n" << pair.second->to_string();
		}
		return stringstream.str();
	}

	template <typename T>
	T get_input(std::string output)
	{
		T value;
		cout << string(50, '\n');
		cout << output;
		cin.clear();
		cin.ignore();
		cin >> value;
		return value;
	}

	void add_matrix(CMtx<TYPE>* mtx)
	{
		string name(get_input<std::string>("Podaj nazwę macierzy: "));
		if (_matrixes.count(name) != 0)
			delete _matrixes.at(name);
		_matrixes.insert_or_assign(name, mtx);
	}

	void add_optional_matrix(COptional<CMtx<TYPE>> optional_mtx)
	{
		if (optional_mtx.is_correct())
		{
			add_matrix(optional_mtx.get_value());
		}
		else
		{
			cout << get_message_for_code(optional_mtx.get_code()) << "\n";
			_getch();
		}
	}

	COptional<int> get_rows()
	{
		COptional<int> optional;
		int* row_cnt = new int(get_input<int>("Podaj ilość wierszy: "));
		if (*row_cnt > 0)
			optional.set_value(row_cnt);
		else
		{
			optional.set_code(CODE_BAD_VALUE);
		}
		return optional;
	}

	COptional<int> get_cols()
	{
		COptional<int> optional;
		int* col_cnt = new int(get_input<int>("Podaj ilość kolumn: "));
		if (*col_cnt > 0)
			optional.set_value(col_cnt);
		else
		{
			optional.set_code(CODE_BAD_VALUE);
		}
		return optional;
	}

	COptional<CMtx<TYPE>> new_matrix_from_console()
	{
		COptional<CMtx<TYPE>> optional;
		COptional<int> optional_rows = get_rows();
		COptional<int> optional_cols = get_cols();

		if (!optional_rows.is_correct())
			optional.set_code(optional_rows.get_code());
		else if (!optional_cols.is_correct())
			optional.set_code(optional_cols.get_code());
		else
		{
			int rows = *optional_rows.get_value();
			int cols = *optional_cols.get_value();
			CMtx<TYPE>* matrix = new CMtx<TYPE>(rows, cols, false);
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++)
				{
					std::stringstream ss;
					ss << "matrix[" << i << "][" << j << "] = ";
					*(*matrix)(i, j) = get_input<TYPE>(ss.str());
				}
			optional.set_value(matrix);
		}
		delete optional_rows.get_value();
		delete optional_cols.get_value();
		return optional;
	}

	COptional<CMtx<TYPE>> new_from_file()
	{
		const std::string file_name(get_input<std::string>("Podaj nazwę pliku: "));
		COptional<CMtx<TYPE>> optional = from_file<TYPE>(file_name);
		return optional;
	}

	COptional<CMtx<TYPE>> new_random_matrix()
	{
		COptional<CMtx<TYPE>> optional;
		COptional<int> optional_rows = get_rows();
		COptional<int> optional_cols = get_cols();
		if (!optional_rows.is_correct())
			optional.set_code(optional_rows.get_code());
		else if (!optional_cols.is_correct())
			optional.set_code(optional_cols.get_code());
		else
			optional.set_value(new CMtx<TYPE>(*optional_rows.get_value(), *optional_cols.get_value(), true));
		delete optional_rows.get_value();
		delete optional_cols.get_value();
		return optional;
	}

	COptional<CMtx<TYPE>> new_matrix()
	{
		switch (get_input<int>("Wybierz opcję : \n1. Losowa\n2. Z pliku\n3. Podana w konsoli\n"))
		{
		case 1:
			return new_random_matrix();
		case 2:
			return new_from_file();
		case 3:
			return new_matrix_from_console();
		default:
			cout << "Niepoprawna wartość\n";
			return new_matrix();
		}
	}


	CMtx<TYPE>* choose_matrix()
	{
		string name(get_input<std::string>(to_string_all_matrixes() + "Podaj nazwę wybranej macierzy: "));
		if (_matrixes.count(name) == 0)
		{
			cout << "\nMacierzy nie znaleziono\n";
			return choose_matrix();
		}
		return _matrixes.at(name);
	};

	template <typename T>
	T choose_number()
	{
		return get_input<int>("Podaj liczbę: ");
	};

	void matrixes_adding()
	{
		add_optional_matrix(*choose_matrix() + *choose_matrix());
	}

	void matrixes_subtracting()
	{
		add_optional_matrix(*choose_matrix() - *choose_matrix());
	}

	void matrixes_multiplication()
	{
		const COptional<CMtx<TYPE>> result(*choose_matrix() * *choose_matrix());
		add_optional_matrix(result);
	}

	void matrix_and_number_multiplication()
	{
		add_matrix(new CMtx<TYPE>(*choose_matrix() * choose_number<TYPE>()));
	}

	void matrix_powering()
	{
		const COptional<CMtx<TYPE>> result(*choose_matrix() ^ choose_number<int>());
		add_optional_matrix(result);
	}

	void matrix_transposition()
	{
		add_matrix(new CMtx<TYPE>(~*choose_matrix()));
	}

	void matrix_det()
	{
		COptional<TYPE> det_optional = det(*choose_matrix());
		if (!det_optional.is_correct())
			cout << get_message_for_code(det_optional.get_code()) << "\n";
		else
			cout << "Wyznacznik wybranej macierzy: " << *det_optional.get_value() << "\n";
		_getch();
		delete det_optional.get_value();
	}

	void vector_dot_product()
	{
		COptional<TYPE> dot_product_optional = choose_matrix()->dot_product(*choose_matrix());
		if (!dot_product_optional.is_correct())
			cout << get_message_for_code(dot_product_optional.get_code()) << "\n";
		else
			cout << "Iloczyn skalarny wektorów: " << *dot_product_optional.get_value() << "\n";
		_getch();
		delete dot_product_optional.get_value();
	}

	void vector_from_matrix_row()
	{
		CMtx<TYPE>* matrix = choose_matrix();
		const int row_index = choose_number<int>();
		const COptional<CMtx<TYPE>> optional_matrix = matrix->get_row(row_index);
		add_optional_matrix(optional_matrix);
	}

	void vector_from_matrix_col()
	{
		CMtx<TYPE>* matrix = choose_matrix();
		const int col_index = choose_number<int>();
		const COptional<CMtx<TYPE>> optional_matrix = matrix->get_column(col_index);
		add_optional_matrix(optional_matrix);
	}

	void matrix_check(void(CProgram::* func)())
	{
		if (!_matrixes.empty())
			(this ->*func)();
		else
		{
			cout << "Brak dostępnej macierzy\n";
			_getch();
		}
	}


	void menu()
	{
		string output = to_string_all_matrixes() +
			"\n\nWybierz opcję : \n1 Nowa macierz\n2 Dodawanie macierzy\n3 Odejmowanie macierzy\n4 Mnożenie macierzy przez macierz\n5 Mnożenie macierzy przez liczbę\n6 Potęgowanie macierzy\n7 Transpozycja macierzy\n8 Wyznacznik macierzy\n9 Iloczyn skalarny (tylko dla macierzy 1naX Xna1\n10 Nowy wektor z wiersza macierzy\n11 Nowy wektor z kolumny macierzy\n12 Wyjście\n Wybór: ";
		switch (get_input<int>(output))
		{
		case 1:
			add_optional_matrix(new_matrix());
			break;
		case 2:
			matrix_check(&CProgram::matrixes_adding);
			break;
		case 3:
			matrix_check(&CProgram::matrixes_subtracting);
			break;
		case 4:
			matrix_check(&CProgram::matrixes_multiplication);
			break;
		case 5:
			matrix_check(&CProgram::matrix_and_number_multiplication);
			break;
		case 6:
			matrix_check(&CProgram::matrix_powering);
			break;
		case 7:
			matrix_check(&CProgram::matrix_transposition);
			break;
		case 8:
			matrix_check(&CProgram::matrix_det);
			break;
		case 9:
			matrix_check(&CProgram::vector_dot_product);
			break;
		case 10:
			matrix_check(&CProgram::vector_from_matrix_row);
			break;
		case 11:
			matrix_check(&CProgram::vector_from_matrix_col);
			break;
		case 12:
			exit(0);
		default:
			std::cout << "\nPodałeś złą liczbę\n";
			_getch();
		}
	}

public:
	void run()
	{
		while (true)
		{
			menu();
		}
	}
};

int main(int argc, char* argv[])
{
	CProgram program;
	program.run();
}
