#include "pch.h"

#include <iostream>
#include <xlocmon>

#include "CppUnitTest.h"
#include "../matrix_default/CMtx.h"
#include "../matrix_default/COptional.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace MyAlgebra;
using namespace MyOptional;

namespace UnitTestMatrixes
{
	TEST_CLASS(UnitTestMatrixes)
	{
	public:
		TEST_METHOD(TestToString)
		{
			CMtx<float> matrix(2, (2), false);
			*matrix(0, 0) = 1.1f;
			*matrix(0, 1) = 2.02f;
			*matrix(1, 0) = 3.003f;
			*matrix(1, 1) = 4.0004f;
			std::string str = matrix.to_string();
			Assert::AreEqual(*new std::string("1.1 2.02 \n3 4 \n"), str);
		}

		TEST_METHOD(TestRandInit)
		{
			int row_cnt = 2;
			int col_cnt = 2;
			CMtx<float> matrix(row_cnt, col_cnt, true);

			for (int i = 0; i < row_cnt; i++)
				for (int j = 0; j < col_cnt; j++)
					for (int l = 0; l < row_cnt; l++)
						for (int k = 0; k < col_cnt; k++)
							if (i != l && j != k)
							{
								std::string string = matrix.to_string();
								Assert::AreNotEqual(*matrix(i, j), *matrix(l, k), 0.000001f,
								                    std::wstring(string.begin(), string.end()).c_str());
							}
		}

		TEST_METHOD(TestDiagonal)
		{
			int row_cnt = 3;
			float value = 6.1255212f;
			CMtx<float> matrix(row_cnt, value);
			for (int i = 0; i < row_cnt; i++)
				for (int j = 0; j < row_cnt; j++)
					if (i == j)
						Assert::AreEqual(value, *matrix(i, j));
					else
						Assert::AreEqual(0.0f, *matrix(i, j));
		}

		TEST_METHOD(TestSettingValues)
		{
			CMtx<float> matrix(2, 2, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(1, 0) = 3;
			*matrix(1, 1) = 4;
			Assert::AreEqual(1.0f, *matrix(0, 0));
			Assert::AreEqual(2.0f, *matrix(0, 1));
			Assert::AreEqual(3.0f, *matrix(1, 0));
			Assert::AreEqual(4.0f, *matrix(1, 1));
		}

		TEST_METHOD(TestAdding)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(2, 2, false);
			*matrix2(0, 0) = 2;
			*matrix2(0, 1) = 4;
			*matrix2(1, 0) = 6;
			*matrix2(1, 1) = 8;

			COptional<CMtx<float>> optional = matrix1 + matrix2;
			Assert::IsTrue(optional.is_correct());
			CMtx<float> matrix3(*optional.get_value());

			Assert::AreEqual(1.0f, *matrix1(0, 0));
			Assert::AreEqual(2.0f, *matrix1(0, 1));
			Assert::AreEqual(3.0f, *matrix1(1, 0));
			Assert::AreEqual(4.0f, *matrix1(1, 1));

			Assert::AreEqual(2.0f, *matrix2(0, 0));
			Assert::AreEqual(4.0f, *matrix2(0, 1));
			Assert::AreEqual(6.0f, *matrix2(1, 0));
			Assert::AreEqual(8.0f, *matrix2(1, 1));

			Assert::AreEqual(3.0f, *matrix3(0, 0));
			Assert::AreEqual(6.0f, *matrix3(0, 1));
			Assert::AreEqual(9.0f, *matrix3(1, 0));
			Assert::AreEqual(12.0f, *matrix3(1, 1));
			delete optional.get_value();
		}

		TEST_METHOD(TestAddingBad)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(1, 2, false);
			*matrix2(0, 0) = 2;
			*matrix2(0, 1) = 4;

			COptional<CMtx<float>> optional(matrix1 + matrix2);
			Assert::AreEqual(CODE_BAD_SIZE, optional.get_code());
		}

		TEST_METHOD(TestAddingBig)
		{
			CMtx<float> matrix1(10, 1, false);
			for (int i = 0; i < 10; i++)
				*matrix1(i, 0) = static_cast<float>(i);

			CMtx<float> matrix2(10, 1, false);
			for (int i = 0; i < 10; i++)
				*matrix2(i, 0) = static_cast<float>(2 * i);

			COptional<CMtx<float>> optional = matrix1 + matrix2;
			Assert::IsTrue(optional.is_correct());
			CMtx<float> result(*optional.get_value());

			for (int i = 0; i < 10; i++)
				Assert::AreEqual(static_cast<float>(3 * i), *result(i, 0));
			delete optional.get_value();
		}

		TEST_METHOD(TestSubstracing)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(2, 2, false);
			*matrix2(0, 0) = 2;
			*matrix2(0, 1) = 4;
			*matrix2(1, 0) = 6;
			*matrix2(1, 1) = 8;

			COptional<CMtx<float>> optional = matrix1 - matrix2;
			Assert::IsTrue(optional.is_correct());
			CMtx<float> matrix3(*optional.get_value());

			Assert::AreEqual(1.0f, *matrix1(0, 0));
			Assert::AreEqual(2.0f, *matrix1(0, 1));
			Assert::AreEqual(3.0f, *matrix1(1, 0));
			Assert::AreEqual(4.0f, *matrix1(1, 1));

			Assert::AreEqual(2.0f, *matrix2(0, 0));
			Assert::AreEqual(4.0f, *matrix2(0, 1));
			Assert::AreEqual(6.0f, *matrix2(1, 0));
			Assert::AreEqual(8.0f, *matrix2(1, 1));

			Assert::AreEqual(-1.0f, *matrix3(0, 0));
			Assert::AreEqual(-2.0f, *matrix3(0, 1));
			Assert::AreEqual(-3.0f, *matrix3(1, 0));
			Assert::AreEqual(-4.0f, *matrix3(1, 1));
			delete optional.get_value();
		}
		
		TEST_METHOD(TestSubstracingBad)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(1, 2, false);
			*matrix2(0, 0) = 2;
			*matrix2(0, 1) = 4;

			COptional<CMtx<float>> optional = matrix1 - matrix2;
			Assert::AreEqual(CODE_BAD_SIZE, optional.get_code());
		}

		TEST_METHOD(TestChangingSign)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = -3;
			*matrix1(1, 1) = 4;

			CMtx<float> result(-matrix1);

			Assert::AreEqual(1.0f, *matrix1(0, 0));
			Assert::AreEqual(2.0f, *matrix1(0, 1));
			Assert::AreEqual(-3.0f, *matrix1(1, 0));
			Assert::AreEqual(4.0f, *matrix1(1, 1));

			Assert::AreEqual(-1.0f, *result(0, 0));
			Assert::AreEqual(-2.0f, *result(0, 1));
			Assert::AreEqual(3.0f, *result(1, 0));
			Assert::AreEqual(-4.0f, *result(1, 1));
		}

		TEST_METHOD(TestComparingFalse)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(2, 2, false);
			*matrix2(0, 0) = 2;
			*matrix2(0, 1) = 4;
			*matrix2(1, 0) = 6;
			*matrix2(1, 1) = 8;

			Assert::IsFalse(matrix1 == matrix2);
		}

		TEST_METHOD(TestComparingTrue)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(2, 2, false);
			*matrix2(0, 0) = 1;
			*matrix2(0, 1) = 2;
			*matrix2(1, 0) = 3;
			*matrix2(1, 1) = 4;

			Assert::IsTrue(matrix1 == matrix2);
		}

		TEST_METHOD(TestAssignCopy)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(2, 2, false);
			*matrix2(0, 0) = 2;
			*matrix2(0, 1) = 4;
			*matrix2(1, 0) = 6;
			*matrix2(1, 1) = 8;

			matrix1.operator=(matrix2);

			Assert::IsTrue(matrix1 == matrix2);
			Assert::IsTrue(~matrix1 == ~matrix2);
		}

		TEST_METHOD(TestAssignDiagonal)
		{
			CMtx<float> matrix(2, 2, false);

			float diagonal = 4.0f;
			matrix.operator=(diagonal);

			Assert::AreEqual(diagonal, *matrix(0, 0));
			Assert::AreEqual(0.0f, *matrix(0, 1));
			Assert::AreEqual(0.0f, *matrix(1, 0));
			Assert::AreEqual(diagonal, *matrix(1, 1));
		}

		TEST_METHOD(TestAssignMove)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(2, 2, false);
			*matrix2(0, 0) = 2;
			*matrix2(0, 1) = 4;
			*matrix2(1, 0) = 6;
			*matrix2(1, 1) = 8;

			CMtx<float> matrix2copy(matrix2);
			CMtx<float> matrix_clear(0, 0, false);

			matrix1.operator=(std::move(matrix2));

			Assert::IsTrue(matrix1 == matrix2copy);
			Assert::IsTrue(~matrix1 == ~matrix2copy);
			Assert::IsTrue(matrix2 == matrix_clear);
		}

		TEST_METHOD(TestNumberMultiplyMatrix)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(1.23f * matrix1);
			Assert::AreEqual(1.23f, *matrix2(0, 0));
			Assert::AreEqual(2.46f, *matrix2(0, 1));
			Assert::AreEqual(3.69f, *matrix2(1, 0));
			Assert::AreEqual(4.92f, *matrix2(1, 1));
		}

		TEST_METHOD(TestMatrixMultiplyNumber)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(matrix1 * 1.23f);
			Assert::AreEqual(1.23f, *matrix2(0, 0));
			Assert::AreEqual(2.46f, *matrix2(0, 1));
			Assert::AreEqual(3.69f, *matrix2(1, 0));
			Assert::AreEqual(4.92f, *matrix2(1, 1));
		}

		TEST_METHOD(TestMultiplyTwoMatrixes)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(2, 2, false);
			*matrix2(0, 0) = 2;
			*matrix2(0, 1) = 4;
			*matrix2(1, 0) = 6;
			*matrix2(1, 1) = 8;

			COptional<CMtx<float>> result_optional = matrix1 * matrix2;

			Assert::AreEqual(CODE_CORRECT, result_optional.get_code());
			CMtx<float>* result = result_optional.get_value();

			Assert::AreEqual(14.0f, *(*result)(0, 0));
			Assert::AreEqual(20.0f, *(*result)(0, 1));
			Assert::AreEqual(30.0f, *(*result)(1, 0));
			Assert::AreEqual(44.0f, *(*result)(1, 1));

			delete result;
		}


		TEST_METHOD(TestMultiplyTwoMatrixesWithDiffrentSizes)
		{
			CMtx<float> matrix1(2, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;

			CMtx<float> matrix2(2, 1, false);
			*matrix2(0, 0) = 2;
			*matrix2(1, 0) = 4;

			COptional<CMtx<float>> result_optional = matrix1 * matrix2;

			Assert::AreEqual(CODE_CORRECT, result_optional.get_code());
			CMtx<float>* result = result_optional.get_value();

			Assert::AreEqual(10.0f, *(*result)(0, 0));
			Assert::AreEqual(22.0f, *(*result)(1, 0));
			delete result;
		}

		TEST_METHOD(TestMultiplyTwoMatrixesWithBadSizes)
		{
			CMtx<float> matrix1(2, 2, false);
			CMtx<float> matrix2(1, 2, false);

			COptional<CMtx<float>> result_optional = matrix1 * matrix2;

			Assert::AreEqual(CODE_BAD_SIZE, result_optional.get_code());
			Assert::IsTrue(nullptr == result_optional.get_value());
		}

		TEST_METHOD(TestMatrixTransposition)
		{
			CMtx<float> matrix1(3, 2, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(1, 0) = 3;
			*matrix1(1, 1) = 4;
			*matrix1(2, 0) = 5;
			*matrix1(2, 1) = 6;

			CMtx<float> matrix2(~matrix1);

			Assert::AreEqual(1.0f, *matrix2(0, 0));
			Assert::AreEqual(3.0f, *matrix2(0, 1));
			Assert::AreEqual(5.0f, *matrix2(0, 2));
			Assert::AreEqual(2.0f, *matrix2(1, 0));
			Assert::AreEqual(4.0f, *matrix2(1, 1));
			Assert::AreEqual(6.0f, *matrix2(1, 2));
		}


		TEST_METHOD(TestMatrixPowerBadValue)
		{
			CMtx<float> matrix1(2, 2, false);

			COptional<CMtx<float>> result = matrix1 ^ -2;

			Assert::AreEqual(CODE_BAD_VALUE, result.get_code());
			Assert::IsTrue(nullptr == result.get_value());
		}

		TEST_METHOD(TestMatrixPowerMinusOne)
		{
			CMtx<float> matrix(2, 2, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(1, 0) = 3;
			*matrix(1, 1) = 4;

			COptional<CMtx<float>> result_optional = matrix ^ -1;

			Assert::AreEqual(CODE_CORRECT, result_optional.get_code());

			CMtx<float> result(*result_optional.get_value());

			Assert::AreEqual(1.0f, *matrix(0, 0));
			Assert::AreEqual(2.0f, *matrix(0, 1));
			Assert::AreEqual(3.0f, *matrix(1, 0));
			Assert::AreEqual(4.0f, *matrix(1, 1));

			Assert::AreEqual(-2.0f, *result(0, 0));
			Assert::AreEqual(1.0f, *result(0, 1));
			Assert::AreEqual(1.5f, *result(1, 0));
			Assert::AreEqual(-0.5f, *result(1, 1));
			delete result_optional.get_value();
		}

		TEST_METHOD(TestMatrixPowerZero)
		{
			CMtx<float> matrix(2, 2, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(1, 0) = 3;
			*matrix(1, 1) = 4;

			COptional<CMtx<float>> result_optional = matrix ^ 0;

			Assert::AreEqual(CODE_CORRECT, result_optional.get_code());

			CMtx<float> result(*result_optional.get_value());

			Assert::AreEqual(1.0f, *matrix(0, 0));
			Assert::AreEqual(2.0f, *matrix(0, 1));
			Assert::AreEqual(3.0f, *matrix(1, 0));
			Assert::AreEqual(4.0f, *matrix(1, 1));

			Assert::AreEqual(1.0f, *result(0, 0));
			Assert::AreEqual(0.0f, *result(0, 1));
			Assert::AreEqual(0.0f, *result(1, 0));
			Assert::AreEqual(1.0f, *result(1, 1));
			delete result_optional.get_value();
		}

		TEST_METHOD(TestMatrixPowerOne)
		{
			CMtx<float> matrix(2, 2, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(1, 0) = 3;
			*matrix(1, 1) = 4;

			COptional<CMtx<float>> result_optional = matrix ^ 1;

			Assert::AreEqual(CODE_CORRECT, result_optional.get_code());

			CMtx<float> result(*result_optional.get_value());

			Assert::AreEqual(1.0f, *matrix(0, 0));
			Assert::AreEqual(2.0f, *matrix(0, 1));
			Assert::AreEqual(3.0f, *matrix(1, 0));
			Assert::AreEqual(4.0f, *matrix(1, 1));

			Assert::AreEqual(1.0f, *result(0, 0));
			Assert::AreEqual(2.0f, *result(0, 1));
			Assert::AreEqual(3.0f, *result(1, 0));
			Assert::AreEqual(4.0f, *result(1, 1));
			delete result_optional.get_value();
		}

		TEST_METHOD(TestMatrixPowerTwo)
		{
			CMtx<float> matrix(2, 2, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(1, 0) = 3;
			*matrix(1, 1) = 4;

			COptional<CMtx<float>> result_optional = matrix ^ 2;

			Assert::AreEqual(CODE_CORRECT, result_optional.get_code());

			CMtx<float> result(*result_optional.get_value());

			Assert::AreEqual(1.0f, *matrix(0, 0));
			Assert::AreEqual(2.0f, *matrix(0, 1));
			Assert::AreEqual(3.0f, *matrix(1, 0));
			Assert::AreEqual(4.0f, *matrix(1, 1));

			Assert::AreEqual(7.0f, *result(0, 0));
			Assert::AreEqual(10.0f, *result(0, 1));
			Assert::AreEqual(15.0f, *result(1, 0));
			Assert::AreEqual(22.0f, *result(1, 1));
			delete result_optional.get_value();
		}

		TEST_METHOD(TestMatrixPowerThree)
		{
			CMtx<float> matrix(2, 2, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(1, 0) = 3;
			*matrix(1, 1) = 4;

			COptional<CMtx<float>> result_optional = matrix ^ 3;

			Assert::AreEqual(CODE_CORRECT, result_optional.get_code());

			CMtx<float>* result = result_optional.get_value();

			Assert::AreEqual(1.0f, *matrix(0, 0));
			Assert::AreEqual(2.0f, *matrix(0, 1));
			Assert::AreEqual(3.0f, *matrix(1, 0));
			Assert::AreEqual(4.0f, *matrix(1, 1));

			Assert::AreEqual(37.0f, *(*result)(0, 0));
			Assert::AreEqual(54.0f, *(*result)(0, 1));
			Assert::AreEqual(81.0f, *(*result)(1, 0));
			Assert::AreEqual(118.0f, *(*result)(1, 1));
			delete result;
		}

		TEST_METHOD(TestScalarProduct)
		{
			CMtx<int> matrix1(1, 4, false);
			*matrix1(0, 0) = 1;
			*matrix1(0, 1) = 2;
			*matrix1(0, 2) = 3;
			*matrix1(0, 3) = 4;

			CMtx<int> matrix2(1, 4, false);
			*matrix2(0, 0) = 2;
			*matrix2(0, 1) = 4;
			*matrix2(0, 2) = 6;
			*matrix2(0, 3) = 8;


			COptional<int> result_optional = matrix1.dot_product(matrix2);

			Assert::AreEqual(CODE_CORRECT, result_optional.get_code());

			int* result = result_optional.get_value();

			Assert::AreEqual(60, *result);
			delete result;
		}

		TEST_METHOD(TestScalarProductColumn)
		{
			CMtx<int> matrix1(4, 1, false);
			*matrix1(0, 0) = 1;
			*matrix1(1, 0) = 2;
			*matrix1(2, 0) = 3;
			*matrix1(3, 0) = 4;

			CMtx<int> matrix2(4, 1, false);
			*matrix2(0, 0) = 2;
			*matrix2(1, 0) = 4;
			*matrix2(2, 0) = 6;
			*matrix2(3, 0) = 8;

			COptional<int> result_optional = matrix1.dot_product(matrix2);
			int* result = result_optional.get_value();

			Assert::AreEqual(CODE_CORRECT, result_optional.get_code());
			Assert::AreEqual(60, *result);
			delete result;
		}

		TEST_METHOD(TestScalarProductBadSize)
		{
			CMtx<int> matrix1(2, 4, false);
			CMtx<int> matrix2(1, 4, false);

			COptional<int> result_optional = matrix1.dot_product(matrix2);

			Assert::AreEqual(CODE_BAD_SIZE, result_optional.get_code());
			Assert::IsTrue(nullptr == result_optional.get_value());
		}

		TEST_METHOD(TestScalarProductBadSize2)
		{
			CMtx<int> matrix1(4, 1, false);
			CMtx<int> matrix2(3, 1, false);

			COptional<int> result_optional = matrix1.dot_product(matrix2);

			Assert::AreEqual(CODE_BAD_SIZE, result_optional.get_code());
			Assert::IsTrue(nullptr == result_optional.get_value());
		}

		TEST_METHOD(TestDet)
		{
			CMtx<int> matrix(2, 2, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(1, 0) = 3;
			*matrix(1, 1) = 4;

			COptional<int> result_optional = det(matrix);
			Assert::AreEqual(CODE_CORRECT, result_optional.get_code());

			int* result = result_optional.get_value();
			Assert::AreEqual(-2, *result);
			delete result;
		}

		TEST_METHOD(TestDetBadSize)
		{
			CMtx<int> matrix(2, 1, false);
			*matrix(0, 0) = 1;
			*matrix(1, 0) = 3;

			COptional<int> result_optional = det(matrix);
			Assert::AreEqual(CODE_BAD_SIZE, result_optional.get_code());
		}

		TEST_METHOD(TestFromFile)
		{
			CMtx<int> matrix(3, 6, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(0, 2) = 3123;
			*matrix(0, 3) = 4;
			*matrix(0, 4) = 5;
			*matrix(0, 5) = 6;
			*matrix(1, 0) = 9;
			*matrix(1, 1) = 11231;
			*matrix(1, 2) = 2;
			*matrix(1, 3) = 3;
			*matrix(1, 4) = 4123;
			*matrix(1, 5) = 6;
			*matrix(2, 0) = 5;
			*matrix(2, 1) = 6;
			*matrix(2, 2) = 723;
			*matrix(2, 3) = 8;
			*matrix(2, 4) = 9;
			*matrix(2, 5) = 9900;

			COptional<CMtx<int>> result = from_file<int>("matrix_int.txt");

			Assert::AreEqual(CODE_CORRECT, result.get_code());
			Assert::IsTrue(matrix == *result.get_value());
			delete result.get_value();
		}

		TEST_METHOD(TestGetDimensions)
		{
			CMtx<int> matrix(1, 2, false);
			Assert::AreEqual(1, matrix.get_row_number());
			Assert::AreEqual(2, matrix.get_col_number());
		}

		TEST_METHOD(TestCreateMatrixBadDimensions)
		{
			CMtx<int> matrix(-1, -2, false);
			Assert::AreEqual(0, matrix.get_row_number());
			Assert::AreEqual(0, matrix.get_col_number());
		}

		TEST_METHOD(TestCreateDiagonalMatrixBadDimensions)
		{
			CMtx<int> matrix(-1, 1);
			Assert::AreEqual(0, matrix.get_row_number());
			Assert::AreEqual(0, matrix.get_col_number());
		}

		TEST_METHOD(TestGetRow)
		{
			CMtx<int> matrix(2, 2, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(1, 0) = 3;
			*matrix(1, 1) = 4;

			COptional<CMtx<int>> optional_row = matrix.get_row(1);
			Assert::IsTrue(optional_row.is_correct());

			CMtx<int> row = *optional_row.get_value();
			Assert::AreEqual(1, row.get_row_number());
			Assert::AreEqual(2, row.get_col_number());
			Assert::AreEqual(3, *row(0, 0));
			Assert::AreEqual(4, *row(0, 1));

			delete optional_row.get_value();
		}

		TEST_METHOD(TestGetRowBad)
		{
			CMtx<int> matrix(2, 2, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(1, 0) = 3;
			*matrix(1, 1) = 4;

			COptional<CMtx<int>> optional_row = matrix.get_row(2);
			Assert::IsFalse(optional_row.is_correct());

			Assert::AreEqual(CODE_BAD_SIZE, optional_row.get_code());
		}
		
		TEST_METHOD(TestGetColumn)
		{
			CMtx<int> matrix(2, 2, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(1, 0) = 3;
			*matrix(1, 1) = 4;

			COptional<CMtx<int>> optional_col = matrix.get_column(1);
			Assert::IsTrue(optional_col.is_correct());

			CMtx<int> col = *optional_col.get_value();
			Assert::AreEqual(2, col.get_row_number());
			Assert::AreEqual(1, col.get_col_number());
			Assert::AreEqual(2, *col(0, 0));
			Assert::AreEqual(4, *col(1, 0));

			delete optional_col.get_value();
		}

		TEST_METHOD(TestGetColumnBad)
		{
			CMtx<int> matrix(2, 2, false);
			*matrix(0, 0) = 1;
			*matrix(0, 1) = 2;
			*matrix(1, 0) = 3;
			*matrix(1, 1) = 4;

			COptional<CMtx<int>> optional_col = matrix.get_column(2);
			Assert::IsFalse(optional_col.is_correct());

			Assert::AreEqual(CODE_BAD_SIZE, optional_col.get_code());
		}
	};
}
