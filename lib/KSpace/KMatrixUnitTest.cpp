
#include <algorithm>
#include <gtest/gtest.h>

#include <KSpace/KMatrix.h>


TEST(KMatrix, Trivial00)
{
    auto a = KMatrix<double>::Ones(3, 6);
    KMatrix<double>::Print(a);

    auto b = KMatrix<double>::Zero(3, 6);
    std::cout << std::endl;
    KMatrix<double>::Print(b);

    double k = 3.;
    auto c = k * KMatrix<double>::Ones(3, 6);
    std::cout << std::endl;
    KMatrix<double>::Print(c);


    auto d = KMatrix<double>::Disc(9, 9, 4, 4, 2);
    std::cout << std::endl;
    KMatrix<double>::Print(d);

    auto e = KMatrix<double>::Line(9, 9, 4, 2, 7);
    std::cout << std::endl;
    KMatrix<double>::Print(e);

    auto h = KMatrix<double>::Rect(9, 9, 2, 2, 5, 5);
    std::cout << std::endl;
    KMatrix<double>::Print(h);

    h[5][5] = 100.;
    h[8][5] = -99.;
    std::cout << std::endl;
    KMatrix<double>::Print(h);

    auto max_vec = KMatrix<double>::Max(h);
    auto min_vec = KMatrix<double>::Min(h);

    auto max_it = std::max_element(max_vec.begin(), max_vec.end());
    auto min_it = std::min_element(min_vec.begin(), min_vec.end());

    std::cout << "max val: " << *max_it << std::endl;
    std::cout << "mix val: " << *min_it << std::endl;
}




