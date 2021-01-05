#include "CVctRef.h"
#include <string>
#include "../CAlgError.h"

namespace RefAlgebra
{
	CVctRef::CVctRef(int size) : size(size)
	{
		column = false;
		vector = new FPTYPE[size];
	}

	CVctRef::CVctRef(const CVctRef& rhs)
	{
		copy(rhs);
	}

	CVctRef::CVctRef(CVctRef&& rhs)
	{
		move(rhs);
	}

	FPTYPE& CVctRef::operator[](int ind) const
	{
		return vector[ind];
	}

	const CVctRef& CVctRef::operator=(const CVctRef& rhs)
	{
		if (&rhs != this)
			copy(rhs);
		return *this;
	}

	const CVctRef& CVctRef::operator=(CVctRef&& rhs)
	{
		if (&rhs != this)
			move(rhs);
		return *this;
	}

	const CVctRef& CVctRef::operator=(FPTYPE val)
	{
		delete vector;
		vector = &val;
		return *this;
	}

	CVctRef CVctRef::map(const CVctRef& rhs, FPTYPE fun(FPTYPE, FPTYPE))
	{
		if (size != rhs.size) throw std::exception("Not allowed");
		CVctRef result(size);
		for (int i = 0; i < size; i++)
		{
			result.vector[i] = fun(vector[i], rhs.vector[i]);
		}
		return result;
	}

	CVctRef CVctRef::operator+(const CVctRef& rhs)
	{
		return map(rhs, [](FPTYPE first, FPTYPE second)
		{
			return first + second;
		});
	}

	CVctRef CVctRef::operator-(const CVctRef& rhs)
	{
		return map(rhs, [](FPTYPE first, FPTYPE second)
		{
			return first - second;
		});
	}

	CVctRef CVctRef::operator*(const CMtxRef& rhs)
	{
		if (rhs.row_cnt != size)
			throw BAD_SIZE_EXCEPTION();
		FPTYPE* result = new FPTYPE[size];
		for (int i = 0; i < size; i++)
		{
			result[i] = 0;
			for (int j = 0; j < rhs.row_cnt; j++)
			{
				result[i] += vector[i] * rhs.row_ptr[i][j];
			}
		}
		return *this;
	}

	CVctRef CVctRef::operator~()
	{
		column = !column;
		return *this;
	}

	bool CVctRef::is_column()
	{
		return column;
	}

	void CVctRef::copy(const CVctRef& rhs)
	{
		delete vector;
		size = rhs.size;
		vector = new FPTYPE[size];
		column = rhs.column;
		for (int i = 0; i < size; i++)
		{
			vector[i] = rhs.vector[i];
		}
	}

	void CVctRef::move(CVctRef& rhs)
	{
		delete vector;
		vector = rhs.vector;
		size = rhs.size;
		column = rhs.column;
		rhs.vector = new FPTYPE[size];
		rhs.column = false;
	}
}
