#pragma once

/*
*	This file contains vector class and its extending functions.
*/

#include <array>
#include <stdexcept>
#include <initializer_list>
#include <cmath>
#include <numeric>
#include <algorithm>

namespace clm
{
	template<int dimensions>
	class Vector
	{
		static_assert(dimensions > 0, "Vector must have at least one dimension!");
		static_assert(dimensions < 5, "Vector class fully supports only 1D-4D vectors.");
	private:
		std::array<float, dimensions> mData;
	public:
		// General ctors

		Vector()
		{
			std::fill(mData.begin(), mData.end(), 0.0f);
		}

		Vector(std::initializer_list<float>&& lst)
		{
			int smallerDim = lst.size() < dimensions ? lst.size() : dimensions;

			for (int i = 0; i < smallerDim; ++i)
			{
				mData[i] = *(lst.begin() + i);
			}
		}

		Vector(float value)
		{
			std::fill(mData.begin(), mData.end(), value);
		}

		// Dimension-specific ctors

		template<size_t dim = dimensions,
			typename std::enable_if_t<dim == 2, int> = 0>
			Vector(float x, float y)
		{
			mData[0] = x;
			mData[1] = y;
		}

		template<size_t dim = dimensions,
			typename std::enable_if_t<dim == 3, int> = 0>
			Vector(float x, float y, float z)
		{
			mData[0] = x;
			mData[1] = y;
			mData[2] = z;
		}

		template<size_t dim = dimensions,
			typename std::enable_if_t<dim == 4, int> = 0>
			Vector(float x, float y, float z, float w)
		{
			mData[0] = x;
			mData[1] = y;
			mData[2] = z;
			mData[3] = w;
		}

		// General setters/getters

		void setComponent(int idx, float value)
		{
			if (idx >= 0 && idx < dimensions)
			{
				mData[idx] = value;
				return;
			}
			throw std::out_of_range("Vector index out of range!");
		}

		float getComponent(int idx) const
		{
			if (idx >= 0 && idx < dimensions)
			{
				return mData[idx];
			}
			throw std::out_of_range("Vector index out of range!");
		}

		int getDimension() const
		{
			return dimensions;
		}

		// Dimension-specific setters/getters

		float x() const
		{
			return mData[0];
		}

		void x(float value)
		{
			mData[0] = value;
		}

		float y() const
		{
			return mData[1];
		}

		void y(float value)
		{
			mData[1] = value;
		}

		template<typename T = float>
		typename std::enable_if_t<std::min(2, dimensions) != dimensions, T>
			z() const
		{
			return mData[2];
		}

		template<typename T = void>
		typename std::enable_if_t<std::min(2, dimensions) != dimensions, T>
			z(float value)
		{
			mData[2] = value;
		}

		template<typename T = float>
		typename std::enable_if_t<std::min(3, dimensions) != dimensions, T>
			w() const
		{
			return mData[3];
		}

		template<typename T = void>
		typename std::enable_if_t<std::min(3, dimensions) != dimensions, T>
			w(float value)
		{
			mData[3] = value;
		}

		// Vector operations

		float length() const
		{
			float sum = 0.0f;

			for (auto c : mData)
			{
				sum += std::pow(c, 2);
			}

			return std::sqrt(sum);
		}

		void normalize()
		{
			float len = length();

			if (len == 0.0f) { return; }

			for (auto& c : mData)
			{
				c /= len;
			}
		}

		float dot(const Vector& rhs) const
		{
			return std::inner_product(mData.begin(), mData.end(), rhs.mData.begin(), 0.0f);
		}

		// Binary cross product of vectors exists only in three and seven dimensional Euclidean space
		template<size_t dim = dimensions,
			typename std::enable_if_t<dim == 3, int> = 0>
			Vector cross(const Vector& rhs) const
		{
			return Vector(
				y() * rhs.z() - z() * rhs.y(),
				z() * rhs.x() - x() * rhs.z(),
				x() * rhs.y() - y() * rhs.x()
				);
		}

		float getAngleWith(const Vector& rhs) const
		{
			float cos = dot(rhs) / (length() * rhs.length());
			cos = cos > 1.0f ? 1.0f : (cos < -1.0f ? -1.0f : cos);
			float val = std::acos(cos);
			return val;
		}

		void negate()
		{
			std::for_each(mData.begin(), mData.end(), [](float& value) { value *= -1.0f; });
		}

		void projectOnto(const Vector& rhs)
		{
			auto tmp = (dot(rhs) / std::pow(rhs.length(), 2)) * rhs;
			mData = tmp.mData;
		}

		// Operators

		float& operator[](size_t idx)
		{
			if (idx >= 0 && idx < dimensions)
			{
				return mData[idx];
			}
			throw std::out_of_range("Vector index out of range!");
		}

		float operator[](size_t idx) const
		{
			return getComponent(idx);
		}

		template<int dim>
		friend std::ostream& operator<<(std::ostream&, const Vector<dim>&);

		template<int dim>
		friend bool operator==(const Vector<dim>&, const Vector<dim>&);

		template<int dim>
		friend bool operator!=(const Vector<dim>&, const Vector<dim>&);

		template<int dim>
		friend Vector<dim> operator*(Vector<dim>, float);

		template<int dim>
		friend Vector<dim> operator*(float, const Vector<dim>&);

		template<int dim>
		friend Vector<dim> operator/(const Vector<dim>&, float);

		Vector& operator*=(float scalar)
		{
			*this = *this * scalar;
			return *this;
		}

		Vector& operator/=(float scalar)
		{
			*this = *this / scalar;
			return *this;
		}

		template<int dim>
		friend Vector<dim> operator+(Vector<dim>, const Vector<dim>&);

		Vector& operator+=(const Vector& rhs)
		{
			*this = *this + rhs;
			return *this;
		}

		template<int dim>
		friend Vector<dim> operator-(Vector<dim>, const Vector<dim>&);

		Vector& operator-=(const Vector& rhs)
		{
			*this = *this - rhs;
			return *this;
		}

	};

	template<int dimensions>
	std::ostream& operator<<(std::ostream& stream, const Vector<dimensions>& vector)
	{
		stream << "(";
		for (size_t i = 0; i < vector.mData.size(); ++i)
		{
			stream << vector.mData[i]; 
			if (i < vector.mData.size() - 1)
			{
				stream << ",";
			}
		}
		stream << ")";

		return stream;
	}

	template<int dimensions>
	bool operator==(const Vector<dimensions>& lhs, const Vector<dimensions>& rhs)
	{
		return std::equal(lhs.mData.begin(), lhs.mData.end(), rhs.mData.begin());
	}

	template<int dimensions>
	bool operator!=(const Vector<dimensions>& lhs, const Vector<dimensions>& rhs)
	{
		return !(lhs == rhs);
	}

	template<int dimensions>
	Vector<dimensions> operator*(Vector<dimensions> vector, float scalar)
	{
		std::for_each(vector.mData.begin(), vector.mData.end(), [&](float& val) { val *= scalar; });
		return vector;
	}

	template<int dimensions>
	Vector<dimensions> operator*(float scalar, const Vector<dimensions>& vector)
	{
		return vector * scalar;
	}

	template<int dimensions>
	Vector<dimensions> operator/(const Vector<dimensions>& vector, float scalar)
	{
		return vector * (1.0f / scalar);
	}

	template<int dimensions>
	Vector<dimensions> operator+(Vector<dimensions> lhs, const Vector<dimensions>& rhs)
	{
		int i = 0;
		std::for_each(lhs.mData.begin(), lhs.mData.end(), [&](float &val) { val += rhs.getComponent(i); ++i; });
		return lhs;
	}

	template<int dimensions>
	Vector<dimensions> operator-(Vector<dimensions> lhs, const Vector<dimensions>& rhs)
	{
		int i = 0;
		std::for_each(lhs.mData.begin(), lhs.mData.end(), [&](float &val) { val -= rhs.getComponent(i); ++i; });
		return lhs;
	}

	// Extending functions
	// -------------------

	template<int dimensions>
	constexpr int getVectorDimensions(const Vector<dimensions>&)
	{
		return dimensions;
	}

	Vector<2> getPerpendicularVector(const Vector<2>& vector)
	{
		return Vector<2>(-vector.y(), vector.x());
	}

	// Resulting vector lies in xy-plane
	Vector<3> getPerpendicularVector(const Vector<3>& vector)
	{
		return Vector<3>(-vector.y(), vector.x(), 0.0f);
	}

	// Resulting vector lies in zw-plane
	Vector<4> getPerpendicularVector(const Vector<4>& vector)
	{
		return Vector<4>(-vector.y(), vector.x(), 0.0f, 0.0f);
	}

	template<int dimensions>
	bool arePointingInSameDirection(const Vector<dimensions>& lhs, const Vector<dimensions>& rhs)
	{
		return lhs.dot(rhs) > 0.0f;
	}

	template<int dimensions>
	bool arePerpendicular(const Vector<dimensions>& lhs, const Vector<dimensions>& rhs)
	{
		return lhs.dot(rhs) == 0.0f;
	}

	template<int dimensions>
	bool areParallel(const Vector<dimensions>& lhs, Vector<dimensions> rhs)
	{
		if (rhs[0] == 0.0f && lhs[0] != 0.0f) { return false; }

		float ratio = lhs[0] / rhs[0];
		rhs *= ratio;

		return lhs == rhs;
	}

	using vec2 = Vector<2>;
	using vec3 = Vector<3>;
	using vec4 = Vector<4>;
}