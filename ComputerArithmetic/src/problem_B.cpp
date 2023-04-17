#include "FPN.h"
#include <iostream>

int main()
{
    class FPN myFPN(2, 3, -1, 1);
    vector_double norm = myFPN.get_all_norm();
    vector_double subnorm = myFPN.get_all_subnorm();

    std::cout << "UFL: " << myFPN.get_UFL() << std::endl;
    std::cout << "OFL: " << myFPN.get_OFL() << std::endl;

    std::cout << "All numbers in F: " << std::endl;
    std::cout <<"0 ";
    for(auto it= norm.begin(); it < norm.end(); it ++)
        std::cout << *it << " " << -*it << " ";
    std::cout << std::endl;

    std::cout << "#F = " << 2 * norm.size() + 1 << std::endl;

    std::cout << "All subnormal numbers of F: " << std::endl;
    for(auto it= subnorm.begin(); it < subnorm.end(); it ++)
        std::cout << *it << " " << -*it << " ";
    std::cout << std::endl;

    return 0;
}