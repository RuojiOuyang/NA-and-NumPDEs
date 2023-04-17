#include <cmath>
#include <iostream>

using namespace std;

#define e 2.718281828
#define PI 3.141592654

int main()
{
    double res =sqrt(2*PI/(1-1/e));
    int i = 2;
    while(1){
        double result = 0;
        double h = 1.0/i;
        result += 1.0*h*pow(e,-0)/2;
        result += 1.0*h*pow(e,-1)/2;
        for (int j = 1; j < i; j++){
            double x = pow(j*h,2);
            result += h*pow(e,-pow(e,-x));
        }
        if (abs(res-result) <= 0.5*pow(10,-6)){
            cout << i;
            break;
        }
    }
    return 0;
}