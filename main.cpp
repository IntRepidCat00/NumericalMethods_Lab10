#include "Integral.h"

int main()
{
//    std::cout << "Enter number of iterations:" << std::endl;
    double eps;
    double a, b;
    std::cout << "Enter eps:" << std::endl;
    std::cin >> eps;
    std::cout << "Enter a border:" << std::endl;
    std::cin >> a;
    std::cout << "Enter b border:" << std::endl;
    std::cin >> b;
    Integral I1(eps, a, b);
    I1.LeftRectanglesMethod(I1.calcNforRc());
    I1.RightRectanglesMethod(I1.calcNforRc());
    I1.MiddleRectanglesMethod(I1.calcNforRc());
    I1.TrapezeMethod(I1.calcNforTr());
    I1.SimpsonMethod(I1.calcNforSm());
    return 0;
}
