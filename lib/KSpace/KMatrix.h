



#ifndef __KMATRIX_H_INCLUDE__
#define __KMATRIX_H_INCLUDE__ 

#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>

template <typename T>

using KMatrixType = std::vector<std::vector<T>>;

template <typename T> 
class KMatrix {
public:
	KMatrix() {};
	KMatrix(KMatrixType<T> &base) { mBase = base; };

	size_t rowSize() { 
		return mBase.size();
	};

	size_t colSize() { 
		if (mBase.empty()) return 0; 
		return mBase.front().size();
	};


	KMatrix operator*(T scalar)
	{
		KMatrix r = *this;

		size_t i = 0;
		for (auto& elem: r.mBase) 
			std::transform(elem.begin(), elem.end(), r.mBase[i++].begin(), 
						   [&] (T &e) { return e * scalar; });
		return r;
	};

	friend KMatrix operator*(T lhs, const KMatrix &rhs) 
	{
		KMatrix<T> r = rhs;

		size_t i = 0;
		for (auto& elem: r.mBase) 
			std::transform(rhs.mBase[i].begin(), rhs.mBase[i].end(), 
						   elem.begin(), [&] (T &e) { return e * lhs; });
		return r;
	}

	friend KMatrix operator*(T lhs, KMatrix &&rhs) 
	{
		for (auto& elem: rhs.mBase) 
			std::transform(elem.begin(), elem.end(), elem.begin(), 
					       [&] (T &e) { return e * lhs; });
		return rhs;
	}

	KMatrix operator+(KMatrix<T>& kmat)
	{
		KMatrix ret = kmat;

		size_t i = 0;
		for (auto& elem: ret.mBase) {
			std::transform(kmat.mBase[i].begin(), kmat.mBase[i].end(), 
						   mBase[i].begin(), elem.begin(), std::plus<T>{});
			i ++;
		}

		return ret;
	}

	std::vector<T> &operator[](int i) {
		return mBase[i];
	}

	static KMatrix Ones(size_t r, size_t c);
	static KMatrix Zero(size_t r, size_t c);
	static KMatrix Disc(size_t r, size_t c, size_t x, size_t y, size_t radius);
	static KMatrix Line(size_t r, size_t c, size_t y, size_t start, size_t end);
	static KMatrix Rect(size_t r, size_t c, size_t left, size_t top, 
			              size_t right, size_t bottom);

	static std::vector<T> Max(KMatrix &mat);
	static std::vector<T> Min(KMatrix &mat);

	static void Print(KMatrix kmat);

private:
	KMatrixType<T> mBase;

};

template <typename T>
KMatrix<T> 
KMatrix<T>::Ones(size_t r, size_t c) 
{
	KMatrix<T> ret;
	ret.mBase = KMatrixType<T>(r, std::vector<T>(c, 1));
	return ret;
}

template <typename T>
KMatrix<T>
KMatrix<T>::Zero(size_t r, size_t c) 
{
	KMatrix<T> ret;
	ret.mBase = KMatrixType<T>(r, std::vector<T>(c, 0));
	return ret;
}

template <typename T>
KMatrix<T>
KMatrix<T>::Disc(size_t r, size_t c, size_t x, size_t y, size_t radius)
{
	assert(x < c && y < r);
	assert(radius < x && x + radius < c);
	assert(radius < y && y + radius < r);

	auto m = Zero(r, c);

	auto x_start = x - radius;
	auto x_end = x + radius;

	auto y_start = y - radius;
	auto y_end = y + radius;
	
	for (size_t i = x_start; i <= x_end; ++ i) 
		for (size_t j = y_start; j <= y_end; ++ j) {
			auto p = i - x;
			auto q = j - y;
			if (p*p + q*q <= radius*radius) m[i][j] = 1;
		}

	return m;
}

template <typename T>
KMatrix<T> 
KMatrix<T>::Line(size_t r, size_t c, size_t y, size_t start, size_t end)
{
	assert(y < r && start < end && end < c);
	auto m = Zero(r, c);

	for (size_t i = start; i <= end; ++ i) m[y][i] = 1;

	return m;
}

template <typename T>
KMatrix<T> 
KMatrix<T>::Rect(size_t r, size_t c, size_t left, size_t top, 
		           size_t right, size_t bottom)
{
	assert(left < right && right < c);
	assert(top < bottom && bottom < r);

	auto m = Zero(r, c);

	for (int i = top; i <= bottom; ++ i) 
		for (int j = left; j <= right; ++ j) 
			m[i][j] = 1;

	return m;
}

template <typename T>
std::vector<T> 
KMatrix<T>::Max(KMatrix<T> &mat)
{
	std::vector<T> ret;

	for (auto &elem : mat.mBase) {
		auto it = std::max_element(elem.begin(), elem.end());
		ret.push_back(*it);
	}
	return ret;
}

template <typename T>
std::vector<T> 
KMatrix<T>::Min(KMatrix<T> &mat)
{
	std::vector<T> ret;

	for (auto &elem : mat.mBase) {
		auto it = std::min_element(elem.begin(), elem.end());
		ret.push_back(*it);
	}
	return ret;
}

template <typename T>
void 
KMatrix<T>::Print(KMatrix<T> kmat)
{
	for (auto &elem: kmat.mBase) {
		std::for_each(elem.begin(), elem.end(), 
					  [](T e) { std::cout << e << " "; });
		std::cout << std::endl;
	}
}

#endif /* ifndef __KMATRIX_H_INCLUDE__ */



