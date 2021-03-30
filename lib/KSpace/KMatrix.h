


#ifndef __KMATRIX_H_INCLUDE__
#define __KMATRIX_H_INCLUDE__

#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <memory>

template<typename T, typename... Args>
T *CreateInstance(Args... args)
{
    return new T(args...);
}

template<typename T, typename... Args>
std::shared_ptr<T> AutoCreateInstance(Args... args)
{
    return std::make_shared<T>(args...);
}

template<typename T>
using KMatrixType = std::vector<std::vector<T>>;

class KBaseMatrix
{
public:
    KBaseMatrix() = default;

    virtual ~KBaseMatrix() = default;
};

template<typename T>
class KMatrix : public KBaseMatrix
{
public:
    KMatrix() = default;

    ~KMatrix() override = default;

    KMatrix(const KMatrix &mat) { mBase = mat.mBase; };

    explicit KMatrix(const KMatrixType<T> &base) { mBase = base; };

    size_t rowSize() const
    {
        return mBase.size();
    };

    size_t colSize() const
    {
        if (mBase.empty()) return 0;
        return mBase.front().size();
    };


    KMatrix operator*(T scalar)
    {
        KMatrix r = *this;

        size_t i = 0;
        for (auto &elem: r.mBase)
            std::transform(elem.begin(), elem.end(), r.mBase[i++].begin(),
                           [&](T &e) { return e * scalar; });
        return r;
    };

    friend KMatrix operator*(T lhs, KMatrix &rhs)
    {
        KMatrix<T> r = rhs;

        size_t i = 0;
        for (auto &elem: r.mBase)
        {
            std::transform(rhs.mBase[i].begin(), rhs.mBase[i].end(),
                           elem.begin(), [=](T &e) { return e * lhs; });
            ++i;
        }
        return r;
    }

    friend KMatrix operator*(T lhs, KMatrix &&rhs)
    {
        for (auto &elem: rhs.mBase)
            std::transform(elem.begin(), elem.end(), elem.begin(),
                           [&](T &e) { return e * lhs; });
        return rhs;
    }

    KMatrix operator+(KMatrix<T> &kmat)
    {
        KMatrix ret = kmat;

        size_t i = 0;
        for (auto &elem: ret.mBase)
        {
            std::transform(kmat.mBase[i].begin(), kmat.mBase[i].end(),
                           mBase[i].begin(), elem.begin(), std::plus<T>{});
            i++;
        }

        return ret;
    }

    KMatrix operator+(KMatrix<T> &&kmat)
    {
        size_t i = 0;
        for (auto &elem: kmat.mBase)
        {
            std::transform(kmat.mBase[i].begin(), kmat.mBase[i].end(),
                           mBase[i].begin(), elem.begin(), std::plus<T>{});
            i++;
        }

        return kmat;
    }

    KMatrix *clone() { return CreateInstance<KMatrix>(*this); }

    std::vector<T> &operator[](int i)
    {
        return mBase[i];
    }

    static KMatrix Ones(size_t r, size_t c);

    static KMatrix Zero(size_t r, size_t c);

    static KMatrix Disc(size_t r, size_t c, size_t x, size_t y, size_t radius);

    static KMatrix Line(size_t r, size_t c, size_t y, size_t start, size_t end);

    static KMatrix Rect(size_t r, size_t c, size_t left, size_t top,
                        size_t right, size_t bottom);

    static void FillDisc(KMatrix<float> &kmat, size_t x, size_t y, size_t radius, T val);

    template<typename F>
    static void Fill(KMatrix<float> &kmat, F lambda);

    static std::vector<T> Max(KMatrix &mat);

    static std::vector<T> Min(KMatrix &mat);

    static void Print(KMatrix kmat);

private:
    KMatrixType<T> mBase;

};

template<typename T>
KMatrix<T>
KMatrix<T>::Ones(size_t r, size_t c)
{
    KMatrix<T> ret;
    ret.mBase = KMatrixType<T>(r, std::vector<T>(c, 1));
    return ret;
}

template<typename T>
KMatrix<T>
KMatrix<T>::Zero(size_t r, size_t c)
{
    KMatrix<T> ret;
    ret.mBase = KMatrixType<T>(r, std::vector<T>(c, 0));
    return ret;
}

template<typename T>
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

    for (size_t i = x_start; i <= x_end; ++i)
        for (size_t j = y_start; j <= y_end; ++j)
        {
            auto p = i - x;
            auto q = j - y;
            if (p * p + q * q <= radius * radius) m[i][j] = 1;
        }

    return m;
}

template<typename T>
KMatrix<T>
KMatrix<T>::Line(size_t r, size_t c, size_t y, size_t start, size_t end)
{
    assert(y < r && start < end && end < c);
    auto m = Zero(r, c);

    for (size_t i = start; i <= end; ++i) m[y][i] = 1;

    return m;
}

template<typename T>
KMatrix<T>
KMatrix<T>::Rect(size_t r, size_t c, size_t left, size_t top,
                 size_t right, size_t bottom)
{
    assert(left < right && right < c);
    assert(top < bottom && bottom < r);

    auto m = Zero(r, c);

    for (size_t i = top; i <= bottom; ++i)
        for (size_t j = left; j <= right; ++j)
            m[i][j] = 1;

    return m;
}

template<typename T>
std::vector<T>
KMatrix<T>::Max(KMatrix<T> &mat)
{
    std::vector<T> ret;

    for (auto &elem : mat.mBase)
    {
        auto it = std::max_element(elem.begin(), elem.end());
        ret.push_back(*it);
    }
    return ret;
}

template<typename T>
std::vector<T>
KMatrix<T>::Min(KMatrix<T> &mat)
{
    std::vector<T> ret;

    for (auto &elem : mat.mBase)
    {
        auto it = std::min_element(elem.begin(), elem.end());
        ret.push_back(*it);
    }
    return ret;
}

template<typename T>
void
KMatrix<T>::Print(KMatrix<T> kmat)
{
    for (auto &elem: kmat.mBase)
    {
        std::for_each(elem.begin(), elem.end(),
                      [](T e) { std::cout << e << " "; });
        std::cout << std::endl;
    }
}

template<typename T>
void
KMatrix<T>::FillDisc(KMatrix<float> &kmat, size_t x, size_t y, size_t radius, T val)
{
    for (size_t i = 0; i < kmat.rowSize(); ++ i)
        for (size_t j = 0; j < kmat.colSize(); ++ j)
            if ((i - x) * (i - x) + (j - y) * (j - y) < radius * radius)
                kmat[i][j] = val;
}

template<typename T>
template<typename F>
void KMatrix<T>::Fill(KMatrix<float> &kmat, F lambda)
{
    for (size_t i = 0; i < kmat.rowSize(); ++ i)
        for (size_t j = 0; j < kmat.colSize(); ++ j)
            lambda(kmat[i][j], i, j);
}

#endif /* ifndef __KMATRIX_H_INCLUDE__ */




