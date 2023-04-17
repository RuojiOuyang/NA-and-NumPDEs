#include "DLS.h"

int main()
{
    vector_double x = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0};
    vector_double y = {2.9, 2.7, 4.8, 5.3, 7.1, 7.6, 7.7, 7.6, 9.4, 9.0, 9.6, 10, 10.2, 9.7, 8.3, 8.4, 9.0, 8.3, 6.6, 6.7, 4.1};
    //vector_double x = {1,2,3,4,5,6,7,8,9,10,11,12};
    //vector_double y = {256,201,159,61,77,40,17,25,103,156,222,345};

    Conditions cond;
    cond.x = x; 
    cond.y = y;
    DLS myDLS;
    myDLS.solve(cond);
    std::cout << myDLS;
    return 0;
}