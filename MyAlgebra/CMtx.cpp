#include "CMtx.h"

#include <cstddef>

#include "CVct.h"

namespace MyAlgebra
{
	CMtx::CMtx(int row_cnt, int col_cnt, bool rand_init)
	{
		this->row_cnt = row_cnt;
		this->col_cnt = col_cnt;
		row_ptr = new FPTYPE * [row_cnt];
		for (int i = 0; i < row_cnt; i++)
			row_ptr[i] = new FPTYPE[col_cnt];
	}

	CMtx::CMtx(int row_cnt, FPTYPE diagonal)
	{
		this->row_cnt = row_cnt;
		this->col_cnt = row_cnt;
		row_ptr = new FPTYPE * [row_cnt];
		for (int i = 0; i < row_cnt; i++)
			row_ptr[i] = new FPTYPE[col_cnt];
	}

	CMtx::CMtx(const CMtx& rhs)
	{
		row_cnt = rhs.row_cnt;
		col_cnt = rhs.col_cnt;
		row_ptr = new FPTYPE * [row_cnt];
		for (int i = 0; i < row_cnt; i++)
		{
			row_ptr[i] = new FPTYPE[col_cnt];
			for(int j = 0; j< col_cnt; j++)
			{
				row_ptr[i][j] = rhs.row_ptr[i][j];
			}
		}
	}

	CMtx::~CMtx()
	{
		if(row_ptr != NULL)
		{
			for(int i = 0; i< row_cnt; i++)
			{
				delete row_ptr[i];
			}
			delete row_ptr;
		}
	}

	const CMtx& CMtx::operator=(const FPTYPE diagonal)
	{
		return *this;
	}

	const CMtx& CMtx::operator=(CMtx&& rhs)
	{
		return *this;
	}

	FPTYPE* CMtx::operator[](int row_ind)
	{
		return row_ptr[row_ind];
	}

	CVct CMtx::operator*(const CVct& rhs) const
	{
		CVct result(10);
		return result;
	}

	CMtx CMtx::operator*(const CMtx& rhs) const
	{
		return *this;
	}

	CMtx CMtx::operator*(FPTYPE multiplier) const
	{
		return *this;
	}

	CMtx CMtx::operator+(const CMtx& rhs) const
	{
		return *this;
	}

	CMtx CMtx::operator-(const CMtx& rhs) const
	{
		return *this;
	}

	CMtx CMtx::operator-() const
	{
		return *this;
	}

	CMtx CMtx::operator~() const
	{
		return *this;
	}

	CMtx CMtx::operator^(int power) const
	{
		return *this;
	}

	bool CMtx::operator==(const CMtx&& rhs) const
	{
		return false;
	}

	void CMtx::display() const
	{
	}

	const CMtx& CMtx::operator=(const CMtx& rhs)
	{
		return *this;
	}
}
