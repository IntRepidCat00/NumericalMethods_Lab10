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
    double calcDeriv1(double x);
    double calcDeriv2(double x);
    double calcDeriv4(double x);
    double RectanglesError(int N);
    double TrapezeError(int N);
    double SimpsonError(int N);
    void SimpsonMethod(int N);
    void LeftRectanglesMethod(int N);
    void RightRectanglesMethod(int N);
    void MiddleRectanglesMethod(int N);
    void TrapezeMethod(int N);
};


#endif //LAB10_INTEGRAL_H
