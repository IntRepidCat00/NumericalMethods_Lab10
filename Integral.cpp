#include "Integral.h"

double Integral::funcAt(double x)
{
    return sqrt(x) * log(1 + cbrt(pow(x, 2)));
}

double Integral::calcDeriv(double x, int order)
{
    const double h = 1.0e-3;

    if (order == 0) return funcAt(x);

    double y1 = calcDeriv(x - h, order - 1);
    double y2 = calcDeriv(x + h, order - 1);

    return (y2 - y1) / (2 * h);
}

double Integral::RectanglesError(int N)
{
    double error;
    double h = (Bbord - Abord) / N;
    double mid = (Bbord - Abord)/2;
    error = (calcDeriv(mid, 2)/24)*(Bbord-Abord)*pow(h, 2);
    return error;
}

double Integral::TrapezeError(int N)
{
    double error;
    double h = (Bbord - Abord) / N;
    double mid = (Bbord - Abord)/2;
    error = -(calcDeriv(mid, 2)/12)*(Bbord-Abord)*pow(h, 2);
    return error;
}

double Integral::SimpsonError(int N)
{
    double error;
    double m;

    double max{calcDeriv(Bbord-Abord, 4)/2};
    if(fabs(max) < fabs(calcDeriv(Abord, 4)))
    {
        max = calcDeriv(Abord, 4);
    }
    if(fabs(max) < fabs(calcDeriv(Bbord, 4)))
    {
        max = calcDeriv(Bbord, 4);
    }

    m=max-(pow(Bbord-Abord, 5))/(180*pow(N, 4));
    error=-(pow(Bbord-Abord, 5)*m)/(180*pow(N, 4));

    return error;
}

double Integral::LeftRectanglesMethod(int N, bool printOrNo)
{
    std::string divideLine(60, '=');
    if(printOrNo)
    {
        std::cout << divideLine << std::endl;
        std::cout << "Left Rectangles Method of Integration" << std::endl;
        std::cout << "N = " << N << std::endl;
        std::cout << "Error = " << RectanglesError(N) << std::endl;
        std::cout << divideLine << std::endl;
    }
    double h{(Bbord-Abord) / N};
    double I{0};

    double Xi{Abord};
    for(int i{0}; i <= N-1; i++)
    {
        double Si{funcAt(Xi)*h};
        I += Si;
        if(printOrNo)
        {
            std::cout << "X[" << i << "] = " << Xi << " | S[" << i << "] = " << Si << std::endl;
        }
        Xi += h;
    }
    if(printOrNo)
    {
        std::cout << "LR | Integral from A = " << Abord << " to B = " << Bbord << " is " << I << std::endl;
    }
    return I;
}

double Integral::RightRectanglesMethod(int N, bool printOrNo)
{
    std::string divideLine(60, '=');
    if(printOrNo)
    {
        std::cout << divideLine << std::endl;
        std::cout << "Right Rectangles Method of Integration" << std::endl;
        std::cout << "N = " << N << std::endl;
        std::cout << "Error = " << RectanglesError(N) << std::endl;
        std::cout << divideLine << std::endl;
    }
    double h{(Bbord-Abord) / N};
    double I{0};

    double Xi{Bbord};
    for(int i{1}; i <= N; i++)
    {
        double Si{funcAt(Xi)*h};
        I += Si;
        if(printOrNo)
        {
            std::cout << "X[" << i << "] = " << Xi << " | S[" << i << "] = " << Si << std::endl;
        }
        Xi -= h;
    }
    if(printOrNo)
    {
        std::cout << divideLine << std::endl;
        std::cout << "RR | Integral from A = " << Abord << " to B = " << Bbord << " is " << I << std::endl;
        std::cout << divideLine << std::endl;
    }

    return I;
}

double Integral::MiddleRectanglesMethod(int N, bool printOrNo)
{
    std::string divideLine(60, '=');
    if(printOrNo)
    {
        std::cout << divideLine << std::endl;
        std::cout << "Middle Rectangles Method of Integration" << std::endl;
        std::cout << "N = " << N << std::endl;
        std::cout << "Error = " << RectanglesError(N) << std::endl;
        std::cout << divideLine << std::endl;
    }
    double h{(Bbord-Abord) / N};
    double I{0};

    double Xi{Abord+(h/2)};
    for(int i{0}; i <= N-1; i++)
    {
        double Si{funcAt(Xi)*h};
        I += Si;
        if(printOrNo)
        {
            std::cout << "X[" << i << "] = " << Xi << " | S[" << i << "] = " << Si << std::endl;
        }
        Xi += h;
    }
    if(printOrNo)
    {
        std::cout << divideLine << std::endl;
        std::cout << "MR | Integral from A = " << Abord << " to B = " << Bbord << " is " << I << std::endl;
        std::cout << divideLine << std::endl;
    }

    return I;
}

double Integral::TrapezeMethod(int N, bool printOrNo)
{
    std::string divideLine(60, '=');
    if(printOrNo)
    {
        std::cout << divideLine << std::endl;
        std::cout << "Trapeze Method of Integration" << std::endl;
        std::cout << "N = " << N << std::endl;
        std::cout << "Error = " << TrapezeError(N) << std::endl;
        std::cout << divideLine << std::endl;
    }
    double h{(Bbord-Abord) / N};
    double I{0};

    double Xi{Abord};
    for(int i{0}; i <= N-1; i++)
    {
        double Si{h*(funcAt(Xi) + funcAt(Xi + h))/2};
        I += Si;
        if(printOrNo)
        {
            std::cout << "X[" << i << "] = " << Xi << " | S[" << i << "] = " << Si << std::endl;
        }
        Xi += h;
    }
    if(printOrNo)
    {
        std::cout << divideLine << std::endl;
        std::cout << "Tr | Integral from A = " << Abord << " to B = " << Bbord << " is " << I << std::endl;
        std::cout << divideLine << std::endl;
    }

    return I;
}

double Integral::SimpsonMethod(int N, bool printOrNo)
{
    std::string divideLine(60, '=');
    if (printOrNo)
    {
        std::cout << divideLine << std::endl;
        std::cout << "Simpson Method of Integration" << std::endl;
        std::cout << "N = " << N << std::endl;
        std::cout << "Error = " << SimpsonError(N) << std::endl;
        std::cout << divideLine << std::endl;
    }
    N *= 2;
    double h{(Bbord - Abord) / N};
    double I{0};

    double Xi{Abord};
    int iter{1};
    for (int i{0}; i <= N - 1; i += 2)
    {

        double Si{(h / 3) * (funcAt(Xi) + 4 * funcAt(Xi + h) + funcAt(Xi + 2 * h))};
        I += Si;
        if (printOrNo)
        {
            std::cout << "X[" << i << "] = " << Xi << " | X[" << i + 1 << "] = " << Xi + h << " | X[" << i + 1 + 1
                      << "] = " << Xi + h + h << " | S[" << iter << "] = " << Si << std::endl;
        }
        Xi += 2 * h;
        iter++;
    }
    if (printOrNo)
    {
        std::cout << divideLine << std::endl;
        std::cout << "Sm | Integral from A = " << Abord << " to B = " << Bbord << " is " << I << std::endl;
        std::cout << divideLine << std::endl;
    }
    return I;
}

int Integral::calcNforRc()
{
    double h{sqrt(eps)};
    double I1{LeftRectanglesMethod(h, false)};
    double I2{I1};


    do
    {
        h /= 2;
        I1 = I2;
        I2 = LeftRectanglesMethod(h, false);
    } while(fabs(I2-I1) > eps);

    return static_cast<int>((Bbord - Abord)/h);
}

int Integral::calcNforTr()
{
    double h{sqrt(eps)};
    double I1{TrapezeMethod(h, false)};
    double I2{I1};


    do
    {
        h /= 2;
        I1 = I2;
        I2 = TrapezeMethod(h, false);
    } while(fabs(I2-I1) > eps);

    return static_cast<int>((Bbord - Abord)/h);
}

int Integral::calcNforSm()
{
    double h{pow(eps, 1.0/4)};
    double I1{SimpsonMethod(h, false)};
    double I2{I1};


    do
    {
        h /= 2;
        I1 = I2;
        I2 = SimpsonMethod(h, false);
    } while(fabs(I2-I1) > eps);

    return static_cast<int>((Bbord - Abord)/h);
}