#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

double f1(double x)
{
    return pow(x,8)-8*pow(x,7)+28*pow(x,6)-56*pow(x,5)+70*pow(x,4)-56*pow(x,3)+28*pow(x,2)-8*x+1;
}

double f2(double x)
{
    return (((((((x - 8)*x + 28)*x - 56)*x + 70)*x - 56)*x + 28)*x - 8)*x + 1;
}

double f3(double x)
{
    return pow(x - 1,8);
}

int main()
{
    vector<double> res1;
    vector<double> res2;
    vector<double> res3;

    for(double i = 0.99; i <= 1.01; i += 0.0002)    
    {
        res1.push_back(f1(i));
        res2.push_back(f2(i));
        res3.push_back(f3(i));
    }

    cout << "Function 1:" << endl;
    for(auto it = res1.begin() ; it < res1.end() ; it ++)
        cout << *it << " ";
    cout << endl;

    cout << "Function 2:" << endl;
    for(auto it = res2.begin() ; it < res2.end() ; it ++)
        cout << *it << " ";
    cout << endl;

    cout << "Function 3:" << endl;
    for(auto it = res3.begin() ; it < res3.end() ; it ++)
        cout << *it << " ";
    cout << endl;
    return 0;
}