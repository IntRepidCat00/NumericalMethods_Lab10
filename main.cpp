#include "Integral.h"

int main()
{
    Integral I1;
    std::cout << "Enter number of iterations:" << std::endl;
    int num;
    std::cin >> num;
    I1.LeftRectanglesMethod(num);
    I1.RightRectanglesMethod(num);
    I1.MiddleRectanglesMethod(num);
    I1.TrapezeMethod(num);
    I1.SimpsonMethod(num);
    return 0;
}
