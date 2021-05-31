
#include <algorithm>
#include <gtest/gtest.h>

#include <KSpace/KMatrix.h>


TEST(KMatrix, Trivial00)
{
    auto a = KMatrix<double>::Ones(3, 6);
//    KMatrix<double>::Print(a);

    auto b = KMatrix<double>::Zero(3, 6);
//    std::cout << std::endl;
//    KMatrix<double>::Print(b);

    double k = 3.;
    auto c = k * KMatrix<double>::Ones(3, 6);
//    std::cout << std::endl;
//    KMatrix<double>::Print(c);


    auto d = KMatrix<double>::Disc(9, 9, 4, 4, 2);
//    std::cout << std::endl;
//    KMatrix<double>::Print(d);

    auto e = KMatrix<double>::Line(9, 9, 4, 2, 7);
//    std::cout << std::endl;
//    KMatrix<double>::Print(e);

    auto h = KMatrix<double>::Rect(9, 9, 2, 2, 5, 5);
//    std::cout << std::endl;
//    KMatrix<double>::Print(h);

    h[5][5] = 100.;
    h[8][5] = -99.;
//    std::cout << std::endl;
//    KMatrix<double>::Print(h);

    auto max_vec = KMatrix<double>::Max(h);
    auto min_vec = KMatrix<double>::Min(h);

    auto max_it = std::max_element(max_vec.begin(), max_vec.end());
    auto min_it = std::min_element(min_vec.begin(), min_vec.end());

//    std::cout << "max val: " << *max_it << std::endl;
//    std::cout << "mix val: " << *min_it << std::endl;
}


TEST(KMatrix, Trivial01)
{
    auto Nx = 16;
    auto Ny = 14;

    /* wedge mask */
    auto r1_mask = KMatrix<double>::Rect(Nx, Ny, 0, 0, Ny-1, Nx/2-1);
    /* part mask */
    auto r2_mask = KMatrix<double>::Rect(Nx, Ny, 0, Nx/2, Ny-1, Nx-1);
    /* wedge */
    auto c1_magnitude = 2700.f;

    /* part */
    auto c2_magnitude = 5800.f;

    /* medium */
    auto C0 = c1_magnitude * r1_mask + c2_magnitude * r2_mask;

    /* flaws definition */
    class PrdFillDisc
    {
    public:
        PrdFillDisc(size_t x, size_t y, size_t radius) : mX(x), mY(y), mRadius(radius), mVal(0) {}
        void setVal(double val) { mVal = val; }
        void operator()(double &e, size_t i, size_t j)
        {
            if ((mX - i) * (mX - i) + (mY - j) * (mY - j) < mRadius * mRadius)
                e = mVal;
        }

    private:
        size_t mX;
        size_t mY;
        double mVal;
        double mRadius;
    };

    auto x_pos = 12;
    auto y_pos = 7;
    auto radius = 2;
    auto flaw_c0_magnitude = 340.f;

    auto fillDisc = PrdFillDisc(x_pos, y_pos, radius);
    fillDisc.setVal(flaw_c0_magnitude);
    KMatrix<double>::Fill(C0, fillDisc);

    KMatrix<double>::Print(C0);
}


