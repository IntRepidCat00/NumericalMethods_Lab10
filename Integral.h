#ifndef LAB10_INTEGRAL_H
#define LAB10_INTEGRAL_H

#include <iostream>
#include <cmath>

class Integral
{
private:
    double eps{10e-3};
    double Abord{0};
    double Bbord{3};
public:
    Integral() = default;
    Integral(double Neps, double NAbord, double NBbord)
        : eps{Neps}, Abord{NAbord}, Bbord{NBbord}
    {}
    double funcAt(double x);
    double calcDeriv(double x, int order);
    double RectanglesError(int N);
    double TrapezeError(int N);
    double SimpsonError(int N);
    double SimpsonMethod(int N, bool printOrNo = true);
    double LeftRectanglesMethod(int N, bool printOrNo = true);
    double RightRectanglesMethod(int N, bool printOrNo = true);
    double MiddleRectanglesMethod(int N, bool printOrNo = true);
    double TrapezeMethod(int N, bool printOrNo = true);
    int calcNforRc();
    int calcNforTr();
    int calcNforSm();
};


#endif //LAB10_INTEGRAL_H
