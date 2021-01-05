#include "CMtxRef.h"

#include <sstream>
#include <cmath>
#include "../CAlgError.h"
#include "CVctRef.h"
#include <xpolymorphic_allocator.h>

namespace RefAlgebra
{
	const FPTYPE CMtxRef::ALG_PRECISION = 0.1f;

	CMtxRef::CMtxRef(): CMtxRef(0, 0, false)
	{
	}

	CMtxRef::CMtxRef(int row_cnt, int col_cnt, bool rand_init)
	{
		this->row_cnt = row_cnt;
		this->col_cnt = col_cnt;
		row_ptr = create_matrix(row_cnt, col_cnt, rand_init);
	}

	CMtxRef::CMtxRef(int row_cnt, FPTYPE diagonal)
	{
		this->row_cnt = row_cnt;
		this->col_cnt = row_cnt;
		row_ptr = create_matrix(row_cnt, col_cnt, diagonal);
	}

	CMtxRef::CMtxRef(const CMtxRef& rhs)
	{
		copy(rhs);
	}

	CMtxRef::CMtxRef(CMtxRef&& rhs)
	{
		move(rhs);
	}

	CMtxRef::~CMtxRef()
	{
		if (row_ptr != NULL)
		{
			for (int i = 0; i < row_cnt; i++)
			{
				delete row_ptr[i];
			}
			delete row_ptr;
		}
	}

	const CMtxRef& CMtxRef::operator=(const FPTYPE diagonal)
	{
		make_diagonal(row_cnt, col_cnt, diagonal, row_ptr);
		return *this;
	}

	FPTYPE** CMtxRef::create_matrix(int row_cnt, int col_cnt)
	{
		FPTYPE** ptr = new FPTYPE*[row_cnt];
		for (int i = 0; i < row_cnt; i++)
		{
			ptr[i] = new FPTYPE[col_cnt];
		}
		return ptr;
	}

	FPTYPE** CMtxRef::create_matrix(int row_cnt, int col_cnt, bool rand_init)
	{
		FPTYPE** ptr = create_matrix(row_cnt, col_cnt);
		if (rand_init)
		{
			for (int i = 0; i < row_cnt; i++)
				for (int j = 0; j < col_cnt; j++)
					ptr[i][j] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
		}
		return ptr;
	}

	void CMtxRef::make_diagonal(int row_cnt, int col_cnt, FPTYPE diagonal, FPTYPE** ptr)
	{
		for (int i = 0; i < row_cnt; i++)
			for (int j = 0; j < col_cnt; j++)
			{
				if (i == j)
					ptr[i][j] = diagonal;
				else
					ptr[i][j] = 0.0f;
			}
	}

	FPTYPE** CMtxRef::create_matrix(int row_cnt, int col_cnt, FPTYPE diagonal)
	{
		FPTYPE** ptr = create_matrix(row_cnt, col_cnt);
		make_diagonal(row_cnt, col_cnt, diagonal, ptr);
		return ptr;
	}

	std::string CMtxRef::to_string()
	{
		std::stringstream stream;
		stream.precision(3);
		for (int i = 0; i < row_cnt; i++)
		{
			for (int j = 0; j < col_cnt; j++)
			{
				stream << matrix[i][j] << " ";
			}
			stream << "\n";
		}
		return stream.str();
	}

	void CMtxRef::change_rows(int first, int second)
	{
		FPTYPE* saved = matrix[first];
		matrix[first] = matrix[second];
		matrix[second] = saved;
	}

	CMtxRef operator*(FPTYPE multiplier, const CMtxRef& rhs)
	{
		CMtxRef result(rhs.row_cnt, rhs.col_cnt, false);
		for (int i = 0; i < rhs.row_cnt; i++)
			for (int j = 0; j < rhs.col_cnt; j++)
				result[i][j] = multiplier * rhs[i][j];
		return result;
	}

	CMtxRef** create_sub_matrixes_for_determinant(const CMtxRef& rhs, const int size)
	{
		CMtxRef** matrixes = new CMtxRef*[size];

		for (int i = 0; i < size; i++)
			matrixes[i] = new CMtxRef(size - 1, size - 1, false);

		for (int i = 0; i < size - 1; i++)
		{
			for (int j = 0; j < size; j++)
			{
				for (int k = 0; k < size; k++)
				{
					if (j < k)
					{
						*(matrixes[k][i][j]) = rhs[i + 1][j];
					}
					else if (j > k)
					{
						*(matrixes[k][i][j - 1]) = rhs[i + 1][j];
					}
				}
			}
		}
		return matrixes;
	}

	float det(const CMtxRef& rhs)
	{
		if (rhs.row_cnt != rhs.col_cnt) throw BAD_SIZE_EXCEPTION();
		const int size = rhs.row_cnt;
		if (size == 0) return 1.0f;

		CMtxRef** matrixes = create_sub_matrixes_for_determinant(rhs, size);

		float result = 0;
		int sign = 1;
		for (int i = 0; i < size; i++)
		{
			result += sign * (rhs[0][i] * det(*matrixes[i]));
			sign *= -1;
		}
		delete matrixes;
		return result;
	}

	void CMtxRef::move(CMtxRef& rhs)
	{
		row_cnt = rhs.row_cnt;
		col_cnt = rhs.col_cnt;
		row_ptr = rhs.row_ptr;
		rhs.row_ptr = create_matrix(row_cnt, col_cnt, false);
	}

	const CMtxRef& CMtxRef::operator=(CMtxRef&& rhs)
	{
		if (&rhs != this)
		{
			delete row_ptr;
			move(rhs);
		}
		return *this;
	}

	FPTYPE*& CMtxRef::operator[](int row_ind) const
	{
		return row_ptr[row_ind];
	}

	CVctRef CMtxRef::operator*(const CVctRef& rhs) const
	{
		return rhs;
	}

	CMtxRef CMtxRef::operator*(const CMtxRef& rhs) const
	{
		if (rhs.row_cnt != col_cnt)
			throw BAD_SIZE_EXCEPTION();
		CMtxRef res(row_cnt, rhs.col_cnt, false);
		for (int i = 0; i < row_cnt; i++)
		{
			for (int j = 0; j < rhs.col_cnt; j++)
			{
				res[i][j] = 0;
				for (int k = 0; k < col_cnt; k++)
					res[i][j] += row_ptr[i][k] * rhs.row_ptr[k][j];
			}
		}
		return res;
	}

	CMtxRef CMtxRef::operator*(FPTYPE multiplier) const
	{
		return multiplier * *this;
	}

	CMtxRef CMtxRef::map(const CMtxRef& rhs, FPTYPE fun(FPTYPE, FPTYPE)) const
	{
		CMtxRef result(col_cnt, row_cnt);
		for (int i = 0; i < row_cnt; i++)
			for (int j = 0; j < col_cnt; j++)
				result[i][j] = fun(row_ptr[i][j], rhs.row_ptr[i][j]);
		return result;
	}

	void CMtxRef::copy(const CMtxRef& rhs)
	{
		row_cnt = rhs.row_cnt;
		col_cnt = rhs.col_cnt;
		row_ptr = new FPTYPE*[row_cnt];
		for (int i = 0; i < row_cnt; i++)
		{
			row_ptr[i] = new FPTYPE[col_cnt];
			for (int j = 0; j < col_cnt; j++)
			{
				row_ptr[i][j] = rhs.row_ptr[i][j];
			}
		}
	}

	CMtxRef CMtxRef::operator+(const CMtxRef& rhs) const
	{
		return map(rhs, [](FPTYPE first, FPTYPE second)
		{
			return first + second;
		});
	}

	CMtxRef CMtxRef::operator-(const CMtxRef& rhs) const
	{
		return map(rhs, [](FPTYPE first, FPTYPE second)
		{
			return first - second;
		});
	}

	CMtxRef CMtxRef::operator-() const
	{
		return -1 * *this;
	}

	CMtxRef CMtxRef::operator~() const
	{
		CMtxRef result(col_cnt, row_cnt);
		for (int i = 0; i < row_cnt; i++)
			for (int j = 0; j < col_cnt; j++)
				result[j][i] = row_ptr[i][j];
		return result;
	}

	CMtxRef CMtxRef::operator^(int power) const
	{
		if (power < -1) throw BAD_VALUE_EXCEPTION();
		if (power == -1)
			return reversed();
		if (power == 0)
		{
			CMtxRef diagonal(row_cnt, 1.0f);
			return diagonal;
		}
		CMtxRef result(*this);
		for (int i = 1; i < power; i++)
		{
			result = std::move(result * *this);
		}
		return result;
	}

	CMtxRef CMtxRef::reversed() const
	{
		if (row_cnt != col_cnt) throw BAD_SIZE_EXCEPTION();
		if (abs(det(*this)) < ALG_PRECISION) throw SINGULAR_MATRIX_EXCEPTION();

		const int size = row_cnt;
		const CMtxRef to_reverse(*this);
		CMtxRef result(size, 1.0f);

		for (int index = 0; index < size; index++)
		{
			const FPTYPE divider = to_reverse[index][index];
			for (int col = 0; col < size; col++)
			{
				to_reverse[index][col] /= divider;
				result[index][col] /= divider;
			}
			for (int row = 0; row < size; row++)
			{
				if (row != index)
				{
					const FPTYPE multiplier = to_reverse[row][index];
					for (int col = 0; col < size; col++)
					{
						to_reverse[row][col] -= to_reverse[index][col] * multiplier;
						result[row][col] -= result[index][col] * multiplier;
					}
				}
			}
		}
		return result;
	}

	bool CMtxRef::operator==(const CMtxRef& rhs) const
	{
		if (row_cnt != rhs.row_cnt || col_cnt != rhs.col_cnt)
			return false;
		for (int i = 0; i < row_cnt; i++)
			for (int j = 0; j < col_cnt; j++)
				if (-ALG_PRECISION < (row_ptr[i][j] - rhs.row_ptr[i][j]) < ALG_PRECISION)
					return false;
		return true;
	}

	void CMtxRef::display() const
	{
	}

	const CMtxRef& CMtxRef::operator=(const CMtxRef& rhs)
	{
		if (&rhs != this)
		{
			delete row_ptr;
			copy(rhs);
		}
		return *this;
	}
}
