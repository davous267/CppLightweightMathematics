#pragma once

/*
*	This file contains matrix class and its extending functions.
*/

#include <array>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <cmath>
#include <sstream>

#include "Vector.h"
#include "Functions.h"

namespace clm
{
	template<int rows, int cols = rows>
	class Matrix
	{
		static_assert(rows > 0 && cols > 0, "Matrix dimensions must be greater than zero.");
		static_assert(rows <= 4 && cols <= 4, "Matrix class fully supports only up to 4x4 matrices.");
	private:
		// Data are stored in row-major order so refering to the array
		// should be done in this way: mData[row][column]
		std::array<std::array<float, cols>, rows> mData;

		void inverse1x1()
		{
			if (determinant1x1() == 0.0f) { return; }
			mData[0][0] = 1.0f / mData[0][0];
		}

		void inverse2x2()
		{
			float det = determinant2x2();
			if (det == 0.0f) { return; }

			// NOTE Below assignment must be done exactly this way (same applies to inv3x3/4x4)
			//		because in "invert" method I call one of the pvt inverse* methods depending on matrix's dimensions
			//		but compiler adds all inverse* methods to the class so if res was just Matrix<2>
			//		I would be given errors during compilation when assigning *this = Matrix<2>
			//		where *this might be even Matrix<1/3/4>
			//		Although this never happens in my code, it might happen and so it causes problems during compilation
			auto res = *this;

			res[0][0] = mData[1][1];
			res[1][1] = mData[0][0];
			res[1][0] = -mData[1][0];
			res[0][1] = -mData[0][1];

			*this = (1.0f / det) * res;
		}

		void inverse3x3()
		{
			float det = determinant3x3();
			if (det == 0.0f) { return; }

			auto res = *this;

			res[0][0] = mData[1][1] * mData[2][2] - mData[1][2] * mData[2][1];
			res[0][1] = mData[0][2] * mData[2][1] - mData[0][1] * mData[2][2];
			res[0][2] = mData[0][1] * mData[1][2] - mData[0][2] * mData[1][1];

			res[1][0] = mData[1][2] * mData[2][0] - mData[1][0] * mData[2][2];
			res[1][1] = mData[0][0] * mData[2][2] - mData[0][2] * mData[2][0];
			res[1][2] = mData[0][2] * mData[1][0] - mData[0][0] * mData[1][2];

			res[2][0] = mData[1][0] * mData[2][1] - mData[1][1] * mData[2][0];
			res[2][1] = mData[0][1] * mData[2][0] - mData[0][0] * mData[2][1];
			res[2][2] = mData[0][0] * mData[1][1] - mData[0][1] * mData[1][0];

			*this = (1.0f / det) * res;
		}

		void inverse4x4()
		{
			float det = determinant4x4();
			if (det == 0.0f) { return; }

			auto res = *this;

			res[0][0] = mData[1][1] * mData[2][2] * mData[3][3] +
						mData[1][2] * mData[2][3] * mData[3][1] +
						mData[1][3] * mData[2][1] * mData[3][2] -
						mData[1][1] * mData[2][3] * mData[3][2] -
						mData[1][2] * mData[2][1] * mData[3][3] -
						mData[1][3] * mData[2][2] * mData[3][1];

			res[0][1] = mData[0][1] * mData[2][3] * mData[3][2] +
						mData[0][2] * mData[2][1] * mData[3][3] +
						mData[0][3] * mData[2][2] * mData[3][1] -
						mData[0][1] * mData[2][2] * mData[3][3] -
						mData[0][2] * mData[2][3] * mData[3][1] -
						mData[0][3] * mData[2][1] * mData[3][2];

			res[0][2] = mData[0][1] * mData[1][2] * mData[3][3] +
						mData[0][2] * mData[1][3] * mData[3][1] +
						mData[0][3] * mData[1][1] * mData[3][2] -
						mData[0][1] * mData[1][3] * mData[3][2] -
						mData[0][2] * mData[1][1] * mData[3][3] -
						mData[0][3] * mData[1][2] * mData[3][1];

			res[0][3] = mData[0][1] * mData[1][3] * mData[2][2] +
						mData[0][2] * mData[1][1] * mData[2][3] +
						mData[0][3] * mData[1][2] * mData[2][1] -
						mData[0][1] * mData[1][2] * mData[2][3] -
						mData[0][2] * mData[1][3] * mData[2][1] -
						mData[0][3] * mData[1][1] * mData[2][2];

			res[1][0] = mData[1][0] * mData[2][3] * mData[3][2] +
						mData[1][2] * mData[2][0] * mData[3][3] +
						mData[1][3] * mData[2][2] * mData[3][0] -
						mData[1][0] * mData[2][2] * mData[3][3] -
						mData[1][2] * mData[2][3] * mData[3][0] -
						mData[1][3] * mData[2][0] * mData[3][2];

			res[1][1] = mData[0][0] * mData[2][2] * mData[3][3] +
						mData[0][2] * mData[2][3] * mData[3][0] +
						mData[0][3] * mData[2][0] * mData[3][2] -
						mData[0][0] * mData[2][3] * mData[3][2] -
						mData[0][2] * mData[2][0] * mData[3][3] -
						mData[0][3] * mData[2][2] * mData[3][0];

			res[1][2] = mData[0][0] * mData[1][3] * mData[3][2] +
						mData[0][2] * mData[1][0] * mData[3][3] +
						mData[0][3] * mData[1][2] * mData[3][0] -
						mData[0][0] * mData[1][2] * mData[3][3] -
						mData[0][2] * mData[1][3] * mData[3][0] -
						mData[0][3] * mData[1][0] * mData[3][2];

			res[1][3] = mData[0][0] * mData[1][2] * mData[2][3] +
						mData[0][2] * mData[1][3] * mData[2][0] +
						mData[0][3] * mData[1][0] * mData[2][2] -
						mData[0][0] * mData[1][3] * mData[2][2] -
						mData[0][2] * mData[1][0] * mData[2][3] -
						mData[0][3] * mData[1][2] * mData[2][0];

			res[2][0] = mData[1][0] * mData[2][1] * mData[3][3] +
						mData[1][1] * mData[2][3] * mData[3][0] +
						mData[1][3] * mData[2][0] * mData[3][1] -
						mData[1][0] * mData[2][3] * mData[3][1] -
						mData[1][1] * mData[2][0] * mData[3][3] -
						mData[1][3] * mData[2][1] * mData[3][0];

			res[2][1] = mData[0][0] * mData[2][3] * mData[3][1] +
						mData[0][1] * mData[2][0] * mData[3][3] +
						mData[0][3] * mData[2][1] * mData[3][0] -
						mData[0][0] * mData[2][1] * mData[3][3] -
						mData[0][1] * mData[2][3] * mData[3][0] -
						mData[0][3] * mData[2][0] * mData[3][1];

			res[2][2] = mData[0][0] * mData[1][1] * mData[3][3] +
						mData[0][1] * mData[1][3] * mData[3][0] +
						mData[0][3] * mData[1][0] * mData[3][1] -
						mData[0][0] * mData[1][3] * mData[3][1] -
						mData[0][1] * mData[1][0] * mData[3][3] -
						mData[0][3] * mData[1][1] * mData[3][0];

			res[2][3] = mData[0][0] * mData[1][3] * mData[2][1] +
						mData[0][1] * mData[1][0] * mData[2][3] +
						mData[0][3] * mData[1][1] * mData[2][0] -
						mData[0][0] * mData[1][1] * mData[2][3] -
						mData[0][1] * mData[1][3] * mData[2][0] -
						mData[0][3] * mData[1][0] * mData[2][1];

			res[3][0] = mData[1][0] * mData[2][2] * mData[3][1] +
						mData[1][1] * mData[2][0] * mData[3][2] +
						mData[1][2] * mData[2][1] * mData[3][0] -
						mData[1][0] * mData[2][1] * mData[3][2] -
						mData[1][1] * mData[2][2] * mData[3][0] -
						mData[1][2] * mData[2][0] * mData[3][1];

			res[3][1] = mData[0][0] * mData[2][1] * mData[3][2] +
						mData[0][1] * mData[2][2] * mData[3][0] +
						mData[0][2] * mData[2][0] * mData[3][1] -
						mData[0][0] * mData[2][2] * mData[3][1] -
						mData[0][1] * mData[2][0] * mData[3][2] -
						mData[0][2] * mData[2][1] * mData[3][0];

			res[3][2] = mData[0][0] * mData[1][2] * mData[3][1] +
						mData[0][1] * mData[1][0] * mData[3][2] +
						mData[0][2] * mData[1][1] * mData[3][0] -
						mData[0][0] * mData[1][1] * mData[3][2] -
						mData[0][1] * mData[1][2] * mData[3][0] -
						mData[0][2] * mData[1][0] * mData[3][1];

			res[3][3] = mData[0][0] * mData[1][1] * mData[2][2] +
						mData[0][1] * mData[1][2] * mData[2][0] +
						mData[0][2] * mData[1][0] * mData[2][1] -
						mData[0][0] * mData[1][2] * mData[2][1] -
						mData[0][1] * mData[1][0] * mData[2][2] -
						mData[0][2] * mData[1][1] * mData[2][0];

			*this = (1.0f / det) * res;
		}

		float determinant1x1() const
		{
			return mData[0][0];
		}

		float determinant2x2() const
		{
			return mData[0][0] * mData[1][1] - 
				   mData[1][0] * mData[0][1];
		}

		float determinant3x3() const
		{
			return mData[0][0] * (mData[1][1] * mData[2][2] - mData[1][2] * mData[2][1]) -
				   mData[0][1] * (mData[1][0] * mData[2][2] - mData[1][2] * mData[2][0]) +
				   mData[0][2] * (mData[1][0] * mData[2][1] - mData[1][1] * mData[2][0]);
		}

		float determinant4x4() const
		{
			float result = 0.0f;
			size_t row = 0;
			size_t col = 0;

			// Computation of determinant based on Laplace expansion
			// |M| = a * |M\{0,0}| - b * |M\{0,1}|  + c * |M\{0,2}|  - d * |M\{0,3}| 
			for (size_t i = 0; i < 4; ++i)
			{
				col = i;
				float scalar = mData[row][col] * (i % 2 == 0 ? 1.0f : -1.0f);

				Matrix<3, 3> matrix;
				for (size_t r = 1; r < 4; ++r)
				{
					for (size_t c = 0; c < 4; ++c)
					{
						if (c == col) { continue; }
						matrix[r - 1][c > col ? c - 1 : c] = mData[r][c];
					}
				}

				result += scalar * matrix.determinant();
			}

			return result;
		}
	public:
		Matrix()
		{
			for (size_t i = 0; i < rows; ++i)
			{
				std::fill(mData[i].begin(), mData[i].end(), 0.0f);
			}
		}

		Matrix(float value, bool fillOnlyMainDiagonal = true)
		{
			for (size_t r = 0; r < rows; ++r)
			{
				for (size_t c = 0; c < cols; ++c)
				{
					mData[r][c] = (r == c || !fillOnlyMainDiagonal) ? value : 0.0f;
				}
			}
		}

		Matrix(const std::array<std::array<float, cols>, rows>& data)
			:mData(data)
		{ }

		// Basic getters/setters

		int getRowCount() const
		{
			return rows;
		}

		int getColumnCount() const
		{
			return cols;
		}

		void setElement(int row, int col, float value)
		{
			if (row >= 0 && col >= 0 && row < rows && col < cols)
			{
				mData[row][col] = value;
				return;
			}
			throw std::out_of_range("Matrix index out of range!");
		}

		float getElement(int row, int col) const
		{
			if (row >= 0 && col >= 0 && row < rows && col < cols)
			{
				return mData[row][col];
			}
			throw std::out_of_range("Matrix index out of range!");
		}

		// Matrix operations

		Matrix& negate()
		{
			std::for_each(mData.begin(), mData.end(), [](auto& arr) { 
				std::for_each(arr.begin(), arr.end(), [](auto& value) { value *= -1.0f; });
			});

			return *this;
		}

		Matrix& multiplyElementWise(const Matrix& rhs)
		{
			for (size_t r = 0; r < rows; ++r)
			{
				for (size_t c = 0; c < cols; ++c)
				{
					mData[r][c] *= rhs[r][c];
				}
			}

			return *this;
		}

		bool isDiagonal() const
		{
			for (size_t r = 0; r < rows; ++r)
			{
				for (size_t c = 0; c < cols; ++c)
				{
					if (r != c && mData[r][c] != 0.0f)
					{
						return false;
					}
				}
			}
			return true;
		}

		float maxElement() const
		{
			float maxValue = mData[0][0];
			for (size_t r = 0; r < rows; ++r)
			{
				maxValue = std::max(maxValue, *std::max_element(mData[r].begin(), mData[r].end()));
			}
			return maxValue;
		}

		float minElement() const
		{
			float minValue = mData[0][0];
			for (size_t r = 0; r < rows; ++r)
			{
				minValue = std::min(minValue, *std::min_element(mData[r].begin(), mData[r].end()));
			}
			return minValue;
		}

		// Applies unary operation with signature float(float) to every matrix element
		template<typename UnaryOperation>
		Matrix& applyOperation(UnaryOperation op)
		{
			std::for_each(mData.begin(), mData.end(), [=](auto& arr) {
				std::for_each(arr.begin(), arr.end(), [=](auto& value) { value = op(value); });
			});

			return *this;
		}

		template<typename T = Matrix&>
		typename std::enable_if_t<cols == rows, T>
			inverse()
		{
			switch (rows)
			{
			case 1:
				inverse1x1();
				break;
			case 2:
				inverse2x2();
				break;
			case 3:
				inverse3x3();
				break;
			case 4:
				inverse4x4();
				break;
			}

			return *this;
		}

		template<typename T = float>
		typename std::enable_if_t<cols == rows, T>
			determinant() const
		{
			switch (rows)
			{
			case 1:
				return determinant1x1();
				break;
			case 2:
				return determinant2x2();
				break;
			case 3:
				return determinant3x3();
				break;
			}

			return determinant4x4();
		}

		// Operators
	     
		// operator[] does not perform bounds checking
		// and its return type is auto because it differs depending on number of [] you use
		// When accesing data with single [], you will receive array of values in given row
		// When accessing data with [][], the result is value at given position
		const auto& operator[](size_t idx) const
		{
			return mData[idx];
		}

		auto& operator[](size_t idx)
		{
			return mData[idx];
		}

		template<int r, int c>
		friend std::ostream& operator<<(std::ostream&, const Matrix<r, c>&);

		template<int r, int c>
		friend Matrix<r, c> operator+(Matrix<r, c>, const Matrix<r, c>&);

		Matrix& operator+=(const Matrix& rhs)
		{
			*this = *this + rhs;
			return *this;
		}

		template<int r, int c>
		friend Matrix<r, c> operator-(Matrix<r, c>, const Matrix<r, c>&);

		Matrix& operator-=(const Matrix& rhs)
		{
			*this = *this - rhs;
			return *this;
		}

		template<int r, int c>
		friend Matrix<r, c> operator*(float, Matrix<r, c>);

		template<int r, int c>
		friend Matrix<r, c> operator*(const Matrix<r, c>&, float);

		Matrix& operator*=(float value)
		{
			*this = *this * value;
			return *this;
		}

		template<int r, int c>
		friend Matrix<r, c> operator/(Matrix<r, c>, float);

		Matrix& operator/=(float value)
		{
			*this = *this / value;
			return *this;
		}

		template<int r, int c>
		friend bool operator==(const Matrix<r, c>&, const Matrix<r, c>&);

		template<int r, int c>
		friend bool operator!=(const Matrix<r, c>&, const Matrix<r, c>&);

		template<int n, int m, int p>
		friend Matrix<n, p> operator*(const Matrix<n, m>&, const Matrix<m, p>&);

		Matrix& operator*=(const Matrix& rhs)
		{
			*this = *this * rhs;
			return *this;
		}

		template<int r, int c>
		friend Vector<r> operator*(const Matrix<r, c>&, const Vector<c>&);
	};

	template<int rows, int cols>
	std::ostream& operator<<(std::ostream& stream, const Matrix<rows, cols>& matrix)
	{
		// Find number with most digits and use its width as arg to std::setw
		size_t widthSize = 1;
		std::stringstream ss;

		for (size_t r = 0; r < matrix.mData.size(); ++r)
		{
			for (size_t c = 0; c < matrix.mData[r].size(); ++c)
			{
				ss << matrix.mData[r][c];
				size_t length = ss.str().size();

				ss.clear();
				ss.str("");

				widthSize = std::max(widthSize, length);
			}
		}

		// Send formatted matrix output to stream
		for (size_t r = 0; r < matrix.mData.size(); ++r)
		{
			for (size_t c = 0; c < matrix.mData[r].size(); ++c)
			{
				stream << (r == 0 && c == 0 ? "[ " : "  ") << 
					   std::setw(widthSize) << matrix.mData[r][c] << ", ";
			}

			if (r == matrix.mData.size() - 1)
			{
				stream << "]";
			}
			stream << std::endl;
		}

		return stream;
	}

	template<int rows, int cols>
	Matrix<rows, cols> operator+(Matrix<rows, cols> lhs, const Matrix<rows, cols>& rhs)
	{
		for (size_t r = 0; r < rows; ++r)
		{
			for (size_t c = 0; c < cols; ++c)
			{
				lhs[r][c] += rhs[r][c];
			}
		}

		return lhs;
	}

	template<int rows, int cols>
	Matrix<rows, cols> operator-(Matrix<rows, cols> lhs, const Matrix<rows, cols>& rhs)
	{
		for (size_t r = 0; r < rows; ++r)
		{
			for (size_t c = 0; c < cols; ++c)
			{
				lhs[r][c] -= rhs[r][c];
			}
		}

		return lhs;
	}

	template<int rows, int cols>
	Matrix<rows, cols> operator*(float scalar, Matrix<rows, cols> matrix)
	{
		std::for_each(matrix.mData.begin(), matrix.mData.end(), [=](auto& arr) {
			std::for_each(arr.begin(), arr.end(), [=](auto& value) { value *= scalar; });
		});
		return matrix;
	}

	template<int rows, int cols>
	Matrix<rows, cols> operator*(const Matrix<rows, cols>& matrix, float scalar)
	{
		return scalar * matrix;
	}

	template<int rows, int cols>
	Matrix<rows, cols> operator/(Matrix<rows, cols> matrix, float scalar)
	{
		std::for_each(matrix.mData.begin(), matrix.mData.end(), [=](auto& arr) {
			std::for_each(arr.begin(), arr.end(), [=](auto& value) { value /= scalar; });
		});
		return matrix;
	}

	template<int rows, int cols>
	bool operator==(const Matrix<rows, cols>& lhs, const Matrix<rows, cols>& rhs)
	{
		return std::equal(lhs.mData.begin(), lhs.mData.end(), rhs.mData.begin());
	}

	template<int rows, int cols>
	bool operator!=(const Matrix<rows, cols>& lhs, const Matrix<rows, cols>& rhs)
	{
		return !(lhs == rhs);
	}

	template<int n, int m, int p>
	Matrix<n, p> operator*(const Matrix<n, m>& lhs, const Matrix<m, p>& rhs)
	{
		Matrix<n, p> result;
		float tmp;

		for (size_t r = 0; r < n; ++r)
		{
			for (size_t c = 0; c < p; ++c)
			{
				tmp = 0.0f;
				for (size_t k = 0; k < m; ++k)
				{
					tmp += lhs[r][k] * rhs[k][c];
				}
				result[r][c] = tmp;
			}
		}

		return result;
	}

	template<int rows, int cols>
	Vector<rows> operator*(const Matrix<rows, cols>& matrix, const Vector<cols>& vector)
	{
		Vector<rows> result;
		float tmp;

		for (size_t r = 0; r < rows; ++r)
		{
			auto& row = matrix[r];

			tmp = 0.0f;
			for (size_t c = 0; c < cols; ++c)
			{
				tmp += row[c] * vector[c];
			}
			result[r] = tmp;
		}

		return result;
	}

	// Extending functions
	// -------------------
	template<int rows, int cols>
	constexpr int getMatrixRowCount(const Matrix<rows, cols>&)
	{
		return rows;
	}

	template<int rows, int cols>
	constexpr int getMatrixColumnCount(const Matrix<rows, cols>&)
	{
		return cols;
	}

	template<int rows, int cols>
	Matrix<cols, rows> transpose(const Matrix<rows, cols>& matrix)
	{
		Matrix<cols, rows> result;

		for (size_t r = 0; r < rows; ++r)
		{
			for (size_t c = 0; c < cols; ++c)
			{
				result[c][r] = matrix[r][c];
			}
		}

		return result;
	}

	template<int side>
	Matrix<side, side> inverse(Matrix<side, side> matrix)
	{
		return matrix.inverse();
	}

	// Predefined & transformation matrices
	// -------------------
	template<int size>
	auto identityMatrix()
	{
		return Matrix<size, size>(1.0f);
	}

	template<int size>
	auto scale(const Vector<size>& scaleAxisValues)
	{
		static_assert(size < 4, "Only up to 3-dimensional vector is supported by this function.");

		Matrix<size + 1, size + 1> result;

		for (size_t r = 0; r < size; ++r)
		{
			for (size_t c = 0; c < size; ++c)
			{
				if (r == c)
				{
					result[r][c] = scaleAxisValues[r];
				}
			}
		}

		result[size][size] = 1.0f;

		return result;
	}

	template<int size>
	auto translate(const Vector<size>& translation)
	{
		static_assert(size < 4, "Only up to 3-dimensional vector is supported by this function.");

		Matrix<size + 1, size + 1> result(1.0f);

		for (size_t r = 0; r < size; ++r)
		{
			for (size_t c = 0; c < size + 1; ++c)
			{
				if (c == size)
				{
					result[r][c] = translation[r];
				}
			}
		}

		return result;
	}

	Matrix<3, 3> rotate2D(float degrees)
	{
		Matrix<3, 3> result(1.0f);
		result[0][0] = cosd(degrees);
		result[0][1] = -sind(degrees);
		result[1][0] = sind(degrees);
		result[1][1] = cosd(degrees);
		return result;
	}

	Matrix<4, 4> rotateX(float degrees)
	{
		Matrix<4, 4> result(1.0f);
		result[1][1] = cosd(degrees);
		result[1][2] = -sind(degrees);
		result[2][1] = sind(degrees);
		result[2][2] = cosd(degrees);
		return result;
	}

	Matrix<4, 4> rotateY(float degrees)
	{
		Matrix<4, 4> result(1.0f);
		result[0][0] = cosd(degrees);
		result[2][0] = -sind(degrees);
		result[0][2] = sind(degrees);
		result[2][2] = cosd(degrees);
		return result;
	}

	Matrix<4, 4> rotateZ(float degrees)
	{
		Matrix<4, 4> result(1.0f);
		result[0][0] = cosd(degrees);
		result[0][1] = -sind(degrees);
		result[1][0] = sind(degrees);
		result[1][1] = cosd(degrees);
		return result;
	}

	using mat2 = Matrix<2>;
	using mat3 = Matrix<3>;
	using mat4 = Matrix<4>;
}