#include "FPN.h"    

double FPN::get_UFL()
{
    return std::pow(beta, L);
}

double FPN::get_OFL()
{
    return std::pow(beta, U)*(beta - std::pow(beta, 1-p));
}

vector_double FPN::get_all_norm()
{
    vector_double res;
    auto it = res.begin();
    // d0
    for(int i = 1 ; i < beta ; i ++)
        res.push_back(i);
    
    // d1 ... dp-1
    for(int i = 1 ; i <= p-1 ; i ++)
    {
        int n = res.size();
        for(int d = 1 ; d < beta ; d ++)
        {
            double t = d / std::pow(beta, i);
            for(int j = 0 ; j < n ; j ++)
                res.push_back(res[j] + t);
        }
    }

    ///exp
    int n = res.size();
    for(int i = 0 ; i < n ; i ++)
    {
        res[i] = res[i] * std::pow(beta, L);
    }
    for(int k = 0 ; k < U - L ; k ++)
        for(int i = 0 ; i < n ; i ++)
            res.push_back(res[n*k + i]*beta);

    return res;
}

vector_double FPN::get_all_subnorm()
{
    vector_double res;
    auto it = res.begin();
    // d0
    res.push_back(0);
    
    // d1 ... dp-1
    for(int i = 1 ; i <= p-1 ; i ++)
    {
        int n = res.size();
        for(int d = 1 ; d < beta ; d ++)
        {
            double t = d / std::pow(beta, i);
            for(int j = 0 ; j < n ; j ++)
                res.push_back(res[j] + t);
        }
    }

    int n = res.size();
    for(int i = 0 ; i < n ; i ++)
        res[i] = res[i] * std::pow(beta, L);

    res.erase(res.begin());

    return res;
}