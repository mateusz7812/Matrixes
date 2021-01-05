#include "CMtx.h"


#include <sstream>
#include <cmath>

#include "CAlgError.h"

namespace MyAlgebra2
{
	const float CMtx::ALG_PRECISION = 0.1f;

	CMtx::CMtx(int row_cnt, int col_cnt, bool rand_init)
	{
		this->row_cnt = row_cnt;
		this->col_cnt = col_cnt;
		row_ptr = create_matrix(row_cnt, col_cnt, rand_init);
	}

	CMtx::CMtx(int row_cnt, float diagonal)
	{
		this->row_cnt = row_cnt;
		this->col_cnt = row_cnt;
		row_ptr = create_matrix(row_cnt, col_cnt, diagonal);
	}

	CMtx::CMtx(const CMtx& rhs)
	{
		copy(rhs);
	}

	CMtx::CMtx(CMtx&& rhs)
	{
		move(rhs);
	}

	CMtx::~CMtx()
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

	const CMtx& CMtx::operator=(const float diagonal)
	{
		make_diagonal(row_cnt, col_cnt, diagonal, row_ptr);
		return *this;
	}

	float** CMtx::create_matrix(int row_cnt, int col_cnt)
	{
		float** ptr = new float*[row_cnt];
		for (int i = 0; i < row_cnt; i++)
		{
			ptr[i] = new float[col_cnt];
		}
		return ptr;
	}

	float** CMtx::create_matrix(int row_cnt, int col_cnt, bool rand_init)
	{
		float** ptr = create_matrix(row_cnt, col_cnt);
		if (rand_init)
		{
			for (int i = 0; i < row_cnt; i++)
				for (int j = 0; j < col_cnt; j++)
					ptr[i][j] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
		}
		return ptr;
	}

	void CMtx::make_diagonal(int row_cnt, int col_cnt, float diagonal, float** ptr)
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

	float** CMtx::create_matrix(int row_cnt, int col_cnt, float diagonal)
	{
		float** ptr = create_matrix(row_cnt, col_cnt);
		make_diagonal(row_cnt, col_cnt, diagonal, ptr);
		return ptr;
	}

	std::string CMtx::to_string()
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

	void CMtx::change_rows(int first, int second)
	{
		float* saved = matrix[first];
		matrix[first] = matrix[second];
		matrix[second] = saved;
	}

	CMtx operator*(float multiplier, const CMtx& rhs)
	{
		CMtx result(rhs.row_cnt, rhs.col_cnt, false);
		for (int i = 0; i < rhs.row_cnt; i++)
			for (int j = 0; j < rhs.col_cnt; j++)
				result[i][j] = multiplier * rhs[i][j];
		return result;
	}

	CMtx** create_sub_matrixes_for_determinant(const CMtx& rhs, const int size)
	{
		CMtx** matrixes = new CMtx*[size];

		for (int i = 0; i < size; i++)
			matrixes[i] = new CMtx(size - 1, size - 1, false);

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

	float det(const CMtx& rhs)
	{
		if (rhs.row_cnt != rhs.col_cnt) throw BAD_SIZE_EXCEPTION();
		const int size = rhs.row_cnt;
		if (size == 0) return 1.0f;

		CMtx** matrixes = create_sub_matrixes_for_determinant(rhs, size);

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

	void CMtx::move(CMtx& rhs)
	{
		row_cnt = rhs.row_cnt;
		col_cnt = rhs.col_cnt;
		row_ptr = rhs.row_ptr;
		rhs.row_ptr = create_matrix(row_cnt, col_cnt, false);
	}

	const CMtx& CMtx::operator=(CMtx&& rhs)
	{
		if (&rhs != this)
		{
			delete row_ptr;
			move(rhs);
		}
		return *this;
	}

	float*& CMtx::operator[](int row_ind) const
	{
		return row_ptr[row_ind];
	}

	CMtx CMtx::operator*(const CMtx& rhs) const
	{
		if (rhs.row_cnt != col_cnt)
			throw BAD_SIZE_EXCEPTION();
		CMtx res(row_cnt, rhs.col_cnt, false);
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

	CMtx CMtx::operator*(float multiplier) const
	{
		return multiplier * *this;
	}

	CMtx CMtx::map(const CMtx& rhs, float fun(float, float)) const
	{
		CMtx result(col_cnt, row_cnt);
		for (int i = 0; i < row_cnt; i++)
			for (int j = 0; j < col_cnt; j++)
				result[i][j] = fun(row_ptr[i][j], rhs.row_ptr[i][j]);
		return result;
	}

	void CMtx::copy(const CMtx& rhs)
	{
		row_cnt = rhs.row_cnt;
		col_cnt = rhs.col_cnt;
		row_ptr = new float*[row_cnt];
		for (int i = 0; i < row_cnt; i++)
		{
			row_ptr[i] = new float[col_cnt];
			for (int j = 0; j < col_cnt; j++)
			{
				row_ptr[i][j] = rhs.row_ptr[i][j];
			}
		}
	}

	CMtx CMtx::operator+(const CMtx& rhs) const
	{
		return map(rhs, [](float first, float second)
		{
			return first + second;
		});
	}

	CMtx CMtx::operator-(const CMtx& rhs) const
	{
		return map(rhs, [](float first, float second)
		{
			return first - second;
		});
	}

	CMtx CMtx::operator-() const
	{
		return -1 * *this;
	}

	CMtx CMtx::operator~() const
	{
		CMtx result(col_cnt, row_cnt);
		for (int i = 0; i < row_cnt; i++)
			for (int j = 0; j < col_cnt; j++)
				result[j][i] = row_ptr[i][j];
		return result;
	}

	CMtx CMtx::operator^(int power) const
	{
		if (power < -1) throw BAD_VALUE_EXCEPTION();
		if (power == -1)
			return reversed();
		if (power == 0)
		{
			CMtx diagonal(row_cnt, 1.0f);
			return diagonal;
		}
		CMtx result(*this);
		for (int i = 1; i < power; i++)
		{
			result = std::move(result * *this);
		}
		return result;
	}

	CMtx CMtx::reversed() const
	{
		if (row_cnt != col_cnt) throw BAD_SIZE_EXCEPTION();
		if (abs(det(*this)) < ALG_PRECISION) throw SINGULAR_MATRIX_EXCEPTION();

		const int size = row_cnt;
		const CMtx to_reverse(*this);
		CMtx result(size, 1.0f);

		for (int index = 0; index < size; index++)
		{
			const float divider = to_reverse[index][index];
			for (int col = 0; col < size; col++)
			{
				to_reverse[index][col] /= divider;
				result[index][col] /= divider;
			}
			for (int row = 0; row < size; row++)
			{
				if (row != index)
				{
					const float multiplier = to_reverse[row][index];
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

	bool CMtx::operator==(const CMtx& rhs) const
	{
		if (row_cnt != rhs.row_cnt || col_cnt != rhs.col_cnt)
			return false;
		for (int i = 0; i < row_cnt; i++)
			for (int j = 0; j < col_cnt; j++)
			{
				if ( abs(row_ptr[i][j] - rhs.row_ptr[i][j]) > ALG_PRECISION)
					return false;
			}
		return true;
	}

	void CMtx::display() const
	{
	}

	const CMtx& CMtx::operator=(const CMtx& rhs)
	{
		if (&rhs != this)
		{
			delete row_ptr;
			copy(rhs);
		}
		return *this;
	}
}
