//
// Created by shaune on 2021/6/4.
//

#ifndef KSPACESOLVER_KINTERP_H
#define KSPACESOLVER_KINTERP_H

template <typename T>
class BaseInterp
{
public:
    int n, mm, jsav, cor, dj;
    const T *xx, *yy;

public:
    BaseInterp(std::vector<T> &x, const T *y, int m)
            : n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y)
    {
        dj = std::min(1, (int) std::pow((T) n, 0.25));
    }

    T interp(T x)
    {
        int jlo = cor ? hunt(x) : locate(x);
        return rawInterp(jlo, x);
    }

    int locate(const T x);

    int hunt(const T x);

    virtual T rawInterp(int jlo, T x) = 0;

};

template <typename T>
int BaseInterp<T>::locate(const T x)
{
    int ju, jm, jl;
    if (n < 2 || mm < 2 || mm > n)
        throw ("locate size error");
    bool ascnd = (xx[n - 1] >= xx[0]);
    jl = 0;
    ju = n - 1;
    while (ju - jl > 1)
    {
        jm = (ju + jl) >> 1;
        if (x >= xx[jm] == ascnd)
            jl = jm;
        else
            ju = jm;
    }
    cor = abs(jl - jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
}

template <typename T>
int BaseInterp<T>::hunt(const T x)
{
    int jl = jsav, jm, ju, inc = 1;
    if (n < 2 || mm < 2 || mm > n)
        throw ("hunt size error");
    bool ascnd = (xx[n - 1] >= xx[0]);

    if (jl < 0 || jl > n - 1)
    {
        jl = 0;
        ju = n - 1;
    } else
    {
        if (x >= xx[jl] == ascnd)
        {
            for (;;)
            {
                ju = jl + inc;
                if (ju >= n - 1)
                {
                    ju = n - 1;
                    break;
                } else if (x < xx[ju] == ascnd)
                    break;
                else
                {
                    jl = ju;
                    inc += inc;
                }
            }
        } else
        {
            ju = jl;
            for (;;)
            {
                jl = jl - inc;
                if (jl <= 0)
                {
                    jl = 0;
                    break;
                } else if (x >= xx[jl] == ascnd)
                    break;
                else
                {
                    ju = jl;
                    inc += inc;
                }
            }
        }
    }
    while (ju - jl > 1)
    {
        jm = (ju + jl) >> 1;
        if (x >= xx[jm] == ascnd)
            jl = jm;
        else
            ju = jm;
    }
    cor = abs(jl - jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
}

template <typename T>
class LinearInterp : public BaseInterp<T>
{
public:
    T rawInterp(int j, T x) override
    {
        if (this->xx[j] == this->xx[j + 1])
            return this->yy[j];
        else
            return this->yy[j] +
                   ((x - this->xx[j]) / (this->xx[j + 1] - this->xx[j])) * (this->yy[j + 1] - this->yy[j]);
    }

public:
    LinearInterp(std::vector<T> &xv, std::vector<T> &yv)
            : BaseInterp<T>(xv, &yv[0], 2) {}
};

#endif //KSPACESOLVER_KINTERP_H
