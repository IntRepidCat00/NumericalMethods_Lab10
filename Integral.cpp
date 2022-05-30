#include "Integral.h"

double Integral::funcAt(double x)
{
    return sqrt(x) * log(1 + cbrt(pow(x, 2)));
}

double Integral::calcDeriv1(double x)
{
    return ((2*pow(x, 1.0/6))/(3*(pow(x, 2/3) + 1))) + ((log(pow(x, 2.0/3)+1))/(2*sqrt(x)));
}

double Integral::calcDeriv2(double x)
{
    return -(4/(9*pow((pow(x, 2.0/3) + 1), 2)*pow(x, 1.0/6)))+(4/(9*(pow(x, 2.0/3) + 1)*pow(x, 5.0/6)))-(log(pow(x, 2.0/3)+1)/(4*pow(x, 3.0/2)));
}

double Integral::calcDeriv4(double x)
{
    return -(32/(27*pow((pow(x, 2.0/3)+1), 4)*pow(x, 5.0/6)))+(50/(81*pow((pow(x, 2.0/3)+1), 2)*pow(x, 13.0/6)))+(100.0/(81*(pow(x, 2.0/3)+1)*pow(x, 17.0/6)))-((15*log(pow(x, 2.0/3)+1))/(16*pow(x, 7.0/2)));
}

double Integral::RectanglesError(int N)
{
    double error;
    double h = (Bbord - Abord) / N;
    double mid = (Bbord - Abord)/2;
    error = (calcDeriv2(mid)/24)*(Bbord-Abord)*pow(h, 2);
    return error;
}

double Integral::TrapezeError(int N)
{
    double error;
    double h = (Bbord - Abord) / N;
    double mid = (Bbord - Abord)/2;
    error = -(calcDeriv2(mid)/12)*(Bbord-Abord)*pow(h, 2);
    return error;
}

double Integral::SimpsonError(int N)
{
    double error;
    double m;

    double max{fabs(calcDeriv4((Bbord-Abord)/2))};
    if(max < fabs(calcDeriv4(Abord)))
    {
        max = fabs(calcDeriv4(Abord));
    }
    if(max < fabs(calcDeriv4(Bbord)))
    {
        max = fabs(calcDeriv4(Bbord));
    }

    m=max-(pow(Bbord-Abord, 5)/(180*pow(N, 4)));
    error=-(pow(Bbord-Abord, 5)*m)/(180*pow(N, 4));

    return error;
}

void Integral::LeftRectanglesMethod(int N)
{
    std::string divideLine(60, '=');
    std::cout << divideLine << std::endl;
    std::cout << "Left Rectangles Method of Integration" << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "Error = " << RectanglesError(N) << std::endl;
    std::cout << divideLine << std::endl;
    double h{(Bbord-Abord) / N};
    double I{0};

    double Xi{Abord};
    for(int i{0}; i <= N-1; i++)
    {
        double Si{funcAt(Xi)*h};
        I += Si;
        std::cout << "X[" << i <<  "] = " << Xi << " | S[" << i << "] = " << Si << std::endl;
        Xi += h;
    }

    std::cout << "LR | Integral from A = " << Abord << " to B = " << Bbord << " is " << I << std::endl;
}

void Integral::RightRectanglesMethod(int N)
{
    std::string divideLine(60, '=');
    std::cout << divideLine << std::endl;
    std::cout << "Right Rectangles Method of Integration" << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "Error = " << RectanglesError(N) << std::endl;
    std::cout << divideLine << std::endl;
    double h{(Bbord-Abord) / N};
    double I{0};

    double Xi{Bbord};
    for(int i{1}; i <= N; i++)
    {
        double Si{funcAt(Xi)*h};
        I += Si;
        std::cout << "X[" << i <<  "] = " << Xi << " | S[" << i << "] = " << Si << std::endl;
        Xi -= h;
    }
    std::cout << divideLine << std::endl;
    std::cout << "RR | Integral from A = " << Abord << " to B = " << Bbord << " is " << I << std::endl;
    std::cout << divideLine << std::endl;
}

void Integral::MiddleRectanglesMethod(int N)
{
    std::string divideLine(60, '=');
    std::cout << divideLine << std::endl;
    std::cout << "Middle Rectangles Method of Integration" << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "Error = " << RectanglesError(N) << std::endl;
    std::cout << divideLine << std::endl;
    double h{(Bbord-Abord) / N};
    double I{0};

    double Xi{Abord+(h/2)};
    for(int i{0}; i <= N-1; i++)
    {
        double Si{funcAt(Xi)*h};
        I += Si;
        std::cout << "X[" << i <<  "] = " << Xi << " | S[" << i << "] = " << Si << std::endl;
        Xi += h;
    }
    std::cout << divideLine << std::endl;
    std::cout << "MR | Integral from A = " << Abord << " to B = " << Bbord << " is " << I << std::endl;
    std::cout << divideLine << std::endl;
}

void Integral::TrapezeMethod(int N)
{
    std::string divideLine(60, '=');
    std::cout << divideLine << std::endl;
    std::cout << "Trapeze Method of Integration" << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "Error = " << TrapezeError(N) << std::endl;
    std::cout << divideLine << std::endl;
    double h{(Bbord-Abord) / N};
    double I{0};

    double Xi{Abord};
    for(int i{0}; i <= N-1; i++)
    {
        double Si{h*(funcAt(Xi) + funcAt(Xi + h))/2};
        I += Si;
        std::cout << "X[" << i <<  "] = " << Xi << " | S[" << i << "] = " << Si << std::endl;
        Xi += h;
    }
    std::cout << divideLine << std::endl;
    std::cout << "Tr | Integral from A = " << Abord << " to B = " << Bbord << " is " << I << std::endl;
    std::cout << divideLine << std::endl;
}

void Integral::SimpsonMethod(int N)
{
    std::string divideLine(60, '=');
    std::cout << divideLine << std::endl;
    std::cout << "Simpson Method of Integration" << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "Error = " << SimpsonError(N) << std::endl;
    std::cout << divideLine << std::endl;
    double h{(Bbord-Abord) / N};
    double I{0};

    double Xi{Abord};
    int iter{1};
    for(int i{0}; i <= N-1; i+=2)
    {

        double Si{(h/3)*(funcAt(Xi)+4* funcAt(Xi+h) + funcAt(Xi + 2*h))};
        I += Si;
        std::cout << "X[" << i <<  "] = " << Xi << " | X[" << i+1 <<  "] = " << Xi+h << " | X[" << i+1+1 <<  "] = " << Xi+h+h << " | S[" << iter << "] = " << Si << std::endl;
        Xi += 2*h;
        iter++;
    }
    std::cout << divideLine << std::endl;
    std::cout << "Sm | Integral from A = " << Abord << " to B = " << Bbord << " is " << I << std::endl;
    std::cout << divideLine << std::endl;
}