#pragma once

#include <cmath>

/*
*	This file contains useful mathematical functions and definitions.
*/

namespace clm
{
	const double PI = 3.141592653589793;

	template<typename T>
	T degToRad(T deg)
	{
		return static_cast<T>(deg * (PI / 180));
	}

	template<typename T>
	T radToDeg(T rad)
	{
		return static_cast<T>(rad * (180 / PI));
	}

	template<typename T>
	typename std::enable_if_t<std::is_arithmetic<T>::value, int>
	floor(T value)
	{
		return static_cast<int>(std::floor(value));
	}

	template<typename T>
	typename std::enable_if_t<std::is_arithmetic<T>::value, int>
	ceil(T value)
	{
		return static_cast<int>(std::ceil(value));
	}

	template<typename T>
	typename std::enable_if_t<std::is_floating_point<T>::value, T>
	fract(T value)
	{
		T wholePart;
		T fractPart = std::modf(value, &wholePart);

		return fractPart;
	}

	float cosd(float degress)
	{
		return std::cos(degToRad(degress));
	}

	float sind(float degress)
	{
		return std::sin(degToRad(degress));
	}

	float tand(float degrees)
	{
		return std::tan(degToRad(degrees));
	}

	template<typename T> 
	int sign(T value) 
	{
		return (T(0) < value) - (value < T(0));
	}

	template<typename T>
	float mix(T first, T second, float alpha)
	{
		return (1.0f - alpha) * first + alpha * second;
	}
}