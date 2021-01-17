#pragma once
#include <sstream>
#include <string>

#include "Optional.h"

namespace UnitTestMatrixes
{
	class UnitTestMatrixes;
}

namespace MyAlgebra
{
	template <typename T>
	class CMtx
	{
	public:
		static const T ALG_PRECISION;

		CMtx<T>(int row_cnt, int col_cnt, bool rand_init = false);
		CMtx<T>(int row_cnt, T diagonal);
		CMtx<T>(const CMtx<T>& rhs);
		CMtx<T>(CMtx<T>&& rhs);
		~CMtx<T>();

		CMtx<T>& operator=(T diagonal);
		CMtx<T>& operator=(const CMtx<T>& rhs);
		CMtx<T>& operator=(CMtx<T>&& rhs);

		T* operator()(int row_ind, int col_ind);

		Optional<CMtx<T>> operator*(const CMtx<T>& rhs) const;
		CMtx<T> operator*(T multiplier) const;

		CMtx<T> operator+(const CMtx<T>& rhs) const;
		CMtx<T> operator-(const CMtx<T>& rhs) const;

		CMtx<T> operator-() const;
		CMtx<T> operator~() const;

		Optional<CMtx<T>> operator^(int power);

		bool operator==(const CMtx<T>& rhs) const;

		Optional<T> dot_product(const CMtx<T>& rhs) const;

		void display() const;

		std::string to_string() const;

		void multiple_row_by(int row, T multiplier);
		void substract_rows_times(int substracted, int substracting, T times);

		friend CMtx<T> operator*(T multiplier, const CMtx<T>& rhs)
		{
			CMtx<T> result(rhs._row_cnt, rhs._col_cnt);
			for (int i = 0; i < rhs._row_cnt; i++)
				for (int j = 0; j < rhs._col_cnt; j++)
					*result(i, j) = multiplier * rhs._row_ptr[i][j];
			return std::move(result);
		}

		friend Optional<T> det(const CMtx<T>& rhs)
		{
			Optional<T> result_optional;
			if (rhs._row_cnt != rhs._col_cnt)
				result_optional.set_code(CODE_BAD_SIZE);
			else if (rhs._row_cnt == 0)
			{
				T* value = new T();
				*value = static_cast<T>(1);
				result_optional.set_value(value);
			}
			else
			{
				CMtx<T>** matrixes = create_sub_matrixes_for_determinant(rhs);

				T* result = new T();
				*result = static_cast<T>(0);
				int sign = 1;
				for (int i = 0; i < rhs._row_cnt && result_optional.get_code() == CODE_CORRECT; i++)
				{
					Optional<T> det_optional = det(*matrixes[i]);
					if (det_optional.get_code() == CODE_CORRECT)
					{
						*result += static_cast<T>(sign) * rhs._row_ptr[0][i] * *det_optional.get_value();
						delete det_optional.get_value();
						sign *= -1;
					}
					else
						result_optional.set_code(det_optional.get_code());
				}
				delete matrixes;
				result_optional.set_value(result);
			}
			return result_optional;
		}

		Optional<CMtx<T>> reversed();

		friend class UnitTestMatrixes::UnitTestMatrixes;

	private:
		T** _row_ptr;
		int _row_cnt;
		int _col_cnt;

		CMtx<T>();
		CMtx<T>(int row_cnt, int col_cnt, T** row_ptr);

		void remove_matrix_ptr();

		void move(CMtx<T>& rhs);
		void copy(const CMtx<T>& rhs);

		friend CMtx<T>** create_sub_matrixes_for_determinant(const CMtx<T>& rhs)
		{
			const int size = rhs._row_cnt;
			CMtx<T>** matrixes = new CMtx<T>*[size];

			for (int i = 0; i < size; i++)
				matrixes[i] = new CMtx<T>(size - 1, size - 1, false);

			for (int i = 0; i < size - 1; i++)
			{
				for (int j = 0; j < size; j++)
				{
					for (int k = 0; k < size; k++)
					{
						if (j < k)
						{
							*(*matrixes[k])(i, j) = rhs._row_ptr[i + 1][j];
						}
						else if (j > k)
						{
							*(*matrixes[k])(i, j - 1) = rhs._row_ptr[i + 1][j];
						}
					}
				}
			}
			return matrixes;
		}
	};

	template <typename T>
	const T CMtx<T>::ALG_PRECISION = static_cast<T>(0.000001);

	template <typename T>
	CMtx<T>::CMtx(int row_cnt, int col_cnt, bool rand_init) : CMtx<T>()
	{
		this->_row_cnt = row_cnt;
		this->_col_cnt = col_cnt;
		_row_ptr = new T*[row_cnt];
		if (rand_init)
		{
			for (int i = 0; i < row_cnt; i++)
			{
				_row_ptr[i] = new T[col_cnt];
				for (int j = 0; j < col_cnt; j++)
				{
					_row_ptr[i][j] = static_cast<T>(rand()) / static_cast<T>(RAND_MAX);
				}
			}
		}
		else
		{
			for (int i = 0; i < row_cnt; i++)
			{
				_row_ptr[i] = new T[col_cnt];
			}
		}
	}

	template <typename T>
	CMtx<T>::CMtx(int row_cnt, T diagonal) : CMtx<T>()
	{
		this->_row_cnt = row_cnt;
		this->_col_cnt = row_cnt;
		_row_ptr = new T*[row_cnt];
		for (int i = 0; i < row_cnt; i++)
		{
			_row_ptr[i] = new T[_col_cnt];
			for (int j = 0; j < _col_cnt; j++)
			{
				if (i == j)
					_row_ptr[i][j] = diagonal;
				else
					_row_ptr[i][j] = 0.0f;
			}
		}
	}

	template <typename T>
	CMtx<T>::CMtx(const CMtx<T>& rhs) : CMtx<T>()
	{
		copy(rhs);
	}

	template <typename T>
	CMtx<T>::CMtx(CMtx<T>&& rhs) : CMtx<T>()
	{
		move(rhs);
	}

	template <typename T>
	CMtx<T>::~CMtx()
	{
		remove_matrix_ptr();
	}

	template <typename T>
	CMtx<T>::CMtx(): _row_ptr(nullptr)
	{
	}

	template <typename T>
	CMtx<T>::CMtx(int row_cnt, int col_cnt, T** row_ptr): _row_ptr(row_ptr), _row_cnt(row_cnt), _col_cnt(col_cnt)
	{
	}

	template <typename T>
	CMtx<T>& CMtx<T>::operator=(const T diagonal)
	{
		for (int i = 0; i < _row_cnt; i++)
			for (int j = 0; j < _col_cnt; j++)
			{
				if (i == j)
					_row_ptr[i][j] = diagonal;
				else
					_row_ptr[i][j] = 0.0f;
			}
		return *this;
	}

	template <typename T>
	CMtx<T>& CMtx<T>::operator=(CMtx<T>&& rhs)
	{
		if (&rhs != this)
		{
			move(rhs);
		}
		return *this;
	}

	template <typename T>
	CMtx<T>& CMtx<T>::operator=(const CMtx<T>& rhs)
	{
		if (&rhs != this)
		{
			copy(rhs);
		}
		return *this;
	}

	template <typename T>
	T* CMtx<T>::operator()(int row_ind, int col_ind)
	{
		return &_row_ptr[row_ind][col_ind];
	}

	template <typename T>
	Optional<CMtx<T>> CMtx<T>::operator*(const CMtx<T>& rhs) const
	{
		if (rhs._row_cnt != _col_cnt)
			return *new Optional<CMtx>(CODE_BAD_SIZE);
		CMtx<T>* res = new CMtx<T>(_row_cnt, rhs._col_cnt, false);
		const CMtx<T> rhsTransposed(~rhs);
		for (int i = 0; i < _row_cnt; i++)
		{
			for (int j = 0; j < rhs._col_cnt; j++)
			{
				*(*res)(i, j) = static_cast<T>(0);
				for (int k = 0; k < _col_cnt; k++)
				{
					*(*res)(i, j) += _row_ptr[i][k] * rhsTransposed._row_ptr[j][k];
				}
			}
		}
		return *new Optional<CMtx<T>>(res);
	}

	template <typename T>
	CMtx<T> CMtx<T>::operator*(T multiplier) const
	{
		return std::move(multiplier * *this);
	}

	template <typename T>
	CMtx<T> CMtx<T>::operator+(const CMtx<T>& rhs) const
	{
		CMtx<T> result(_row_cnt, _col_cnt);
		for (int i = 0; i < _row_cnt; i++)
		{
			for (int j = 0; j < _col_cnt; j++)
			{
				*result(i, j) = _row_ptr[i][j] + rhs._row_ptr[i][j];
			}
		}
		return std::move(result);
	}

	template <typename T>
	CMtx<T> CMtx<T>::operator-(const CMtx<T>& rhs) const
	{
		CMtx<T> result(_row_cnt, _col_cnt);
		for (int i = 0; i < _row_cnt; i++)
		{
			T* row = _row_ptr[i];
			T* rhs_row = rhs._row_ptr[i];
			for (int j = 0; j < _col_cnt; j++)
			{
				*result(i, j) = row[j] - rhs_row[j];
			}
		}
		return std::move(result);
	}

	template <typename T>
	CMtx<T> CMtx<T>::operator-() const
	{
		return std::move(static_cast<T>(-1) * *this);
	}

	template <typename T>
	CMtx<T> CMtx<T>::operator~() const
	{
		T** transposed = new T*[_col_cnt];
		for (int i = 0; i < _col_cnt; i++)
		{
			transposed[i] = new T[_row_cnt];
			for (int j = 0; j < _row_cnt; j++)
			{
				transposed[i][j] = _row_ptr[j][i];
			}
		}

		return std::move(CMtx<T>(_col_cnt, _row_cnt, transposed));
	}

	template <typename T>
	Optional<CMtx<T>> CMtx<T>::operator^(int power)
	{
		Optional<CMtx<T>> result_optional;
		if (power < -1)
			result_optional.set_code(CODE_BAD_VALUE);
		else if (power == -1)
			result_optional = reversed();
		else if (power == 0)
		{
			CMtx<T>* diagonal = new CMtx<T>(_row_cnt, static_cast<T>(1));
			result_optional.set_value(diagonal);
		}
		else
		{
			CMtx<T>* result = new CMtx<T>(*this);
			for (int i = 1; i < power && result_optional.get_code() == CODE_CORRECT; i++)
			{
				Optional<CMtx<T>> multiplication_optional = *result * *this;
				if (multiplication_optional.get_code() != CODE_CORRECT)
					result_optional.set_code(multiplication_optional.get_code());
				*result = *multiplication_optional.get_value();
			}
			result_optional.set_value(result);
		}
		return result_optional;
	}

	template <typename T>
	bool CMtx<T>::operator==(const CMtx<T>& rhs) const
	{
		if (_row_cnt != rhs._row_cnt || _col_cnt != rhs._col_cnt)
			return false;
		for (int i = 0; i < _row_cnt; i++)
			for (int j = 0; j < _col_cnt; j++)
			{
				if (abs(_row_ptr[i][j] - rhs._row_ptr[i][j]) > ALG_PRECISION)
					return false;
			}
		return true;
	}

	template <typename T>
	Optional<T> CMtx<T>::dot_product(const CMtx<T>& rhs) const
	{
		Optional<T> result_optional;
		if (_col_cnt != 1 && _row_cnt != 1) result_optional.set_code(CODE_BAD_SIZE);
		else if (_col_cnt != rhs._col_cnt || _row_cnt != rhs._row_cnt) result_optional.set_code(CODE_BAD_SIZE);
		else
		{
			T* result = new T();
			*result = static_cast<T>(0);
			for (int i = 0; i < _row_cnt; i++)
				for (int j = 0; j < _col_cnt; j++)
					*result += _row_ptr[i][j] * rhs._row_ptr[i][j];
			result_optional.set_value(result);
		}
		return result_optional;
	}

	template <typename T>
	std::string CMtx<T>::to_string() const
	{
		std::stringstream stream;
		stream.precision(3);
		for (int i = 0; i < _row_cnt; i++)
		{
			for (int j = 0; j < _col_cnt; j++)
			{
				stream << _row_ptr[i][j] << " ";
			}
			stream << "\n";
		}
		return std::move(stream.str());
	}

	template <typename T>
	void CMtx<T>::multiple_row_by(int row, T multiplier)
	{
		for (int col = 0; col < _col_cnt; col++)
		{
			_row_ptr[row][col] *= multiplier;
		}
	}

	template <typename T>
	void CMtx<T>::substract_rows_times(int substracted, int substracting, T times)
	{
		for (int col = 0; col < _col_cnt; col++)
		{
			_row_ptr[substracted][col] -= _row_ptr[substracting][col] * times;
		}
	}

	template <typename T>
	void CMtx<T>::display() const
	{
		std::cout << to_string();
	}

	template <typename T>
	void CMtx<T>::copy(const CMtx<T>& rhs)
	{
		remove_matrix_ptr();
		_row_cnt = rhs._row_cnt;
		_col_cnt = rhs._col_cnt;
		_row_ptr = new T*[_row_cnt];
		for (int i = 0; i < _row_cnt; i++)
		{
			_row_ptr[i] = new T[_col_cnt];
			T* rhs_row = rhs._row_ptr[i];

			for (int j = 0; j < _col_cnt; j++)
			{
				_row_ptr[i][j] = rhs_row[j];
			}
		}
	}

	template <typename T>
	void CMtx<T>::move(CMtx<T>& rhs)
	{
		remove_matrix_ptr();
		_row_cnt = rhs._row_cnt;
		_col_cnt = rhs._col_cnt;
		_row_ptr = rhs._row_ptr;

		rhs._row_cnt = 0;
		rhs._col_cnt = 0;
		rhs._row_ptr = nullptr;
	}

	template <typename T>
	void CMtx<T>::remove_matrix_ptr()
	{
		if (_row_ptr != nullptr)
		{
			for (int i = 0; i < _row_cnt; i++)
			{
				delete _row_ptr[i];
			}
			delete _row_ptr;
			_row_ptr = nullptr;
		}
	}

	template <typename T>
	Optional<CMtx<T>> CMtx<T>::reversed()
	{
		Optional<CMtx<T>> result_optional;
		Optional<T> det_optional = det(*this);
		if (_row_cnt != _col_cnt)
			result_optional.set_code(CODE_BAD_SIZE);
		else if (det_optional.get_code() != CODE_CORRECT)
			result_optional.set_code(det_optional.get_code());
		else if (abs(*det_optional.get_value()) < ALG_PRECISION)
			result_optional.set_code(CODE_SINGULAR_MATRIX);
		else
		{
			CMtx<T> to_reverse(*this);
			CMtx<T>* result = new CMtx<T>(_row_cnt, 1.0f);

			for (int index = 0; index < _row_cnt; index++)
			{
				T multiplier = static_cast<T>(1) / *to_reverse(index, index);
				result->multiple_row_by(index, multiplier);
				to_reverse.multiple_row_by(index, multiplier);
				for (int row = 0; row < _row_cnt; row++)
				{
					if (row != index)
					{
						T times = *to_reverse(row, index);
						to_reverse.substract_rows_times(row, index, times);
						result->substract_rows_times(row, index, times);
					}
				}
			}
			result_optional.set_value(result);
		}
		return result_optional;
	}
}
