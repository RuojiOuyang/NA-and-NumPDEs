#include "test_helper.h"
#include <fstream>
#include <chrono>

using time_point = std::chrono::time_point<std::chrono::high_resolution_clock>;

double sec_duration(time_point &a, time_point &b)
{
    auto d = std::chrono::duration_cast<std::chrono::duration<double>>(b-a);
    // std::cout << d.count()<<"s" << std::endl;
    return d.count();
}

struct combination
{
    interpolation_type inter;
    restriction_type res;
    string name;
}op[4] = {
    {linear, injection, "(linear + injection)"},
    {linear, full_weighting, "(linear + full-weighting)"},
    {quadratic, injection, "(quadratic + injection)"},
    {quadratic, full_weighting, "(quadratic + full-weighting)"}
};

int helper::output(string _title, string _filename, string _type, int _dim, VectorXd _x)
{
    string save_path = string("../plot/data/") + mark + _filename + (".csv");
    _filename = mark+_filename;
    std::ofstream outfile;
    outfile.open(save_path);
    outfile << _title << std::endl << _filename << std::endl << _type << ", " << _dim << std::endl;
    for(int i = 0 ; i < _x.size() - 1; i ++)
        outfile << _x(i) << ", ";
    outfile << _x(_x.size()-1);
    return 0;
}

int helper1d::generate_json()
{
    json_path.clear();
    string prefix = string("../json/") + mark + string("_");
    for(int N = 32 ; N <= 256 ; N *= 2)
    {
        string path = prefix + std::to_string(N) + string(".json");
        json_path.push_back(path);
        class genJson<1> myJson(N);
        myJson.set_fval(f);
        myJson.set_BV(0, 0, u(0));
        if(type == 0)
            myJson.set_BV(1, 0, u(1));
        else
            myJson.set_BV(1, 1, neumann);
        myJson.write_json(path);
    }
    return 0;
}

double helper1d::get_abs_err(const VectorXd& _x, int _norm)
{
    VectorXd e = get_abs_err(_x);
    if(_norm == 0)
        return e.lpNorm<Eigen::Infinity>();
    if(_norm == 1)
        return e.lpNorm<1>();
    if(_norm == 2)
        return e.lpNorm<2>();
    return D_MAX;
}

VectorXd helper1d::get_abs_err(const VectorXd& _x)
{
    int N = _x.size() - 1;
    double h = 1/double(N);
    VectorXd e(N+1);
    for(int i = 0 ; i <= N ; i ++)
        e(i) = std::abs(u(i*h)-_x(i));
    return e;
}

// int helper1d::test_all_type(Multigrid<1> &_mg)
// {
//     VectorXd x;
//     // linear + injection
//     _mg.set_interpolation(linear);
//     _mg.set_restriction(injection)
//     x = _mg.VC_solver(1,1);
//     output("$N=256$ V-cycle (injection + linear)", "VC_injection_linear", "err", 1, get_abs_err(x));
//     x = _mg.FMG_solver(20, 20);
//     output("$N=256$ V-cycle (injection + linear)", "VC_injection_linear", "err", 1, get_abs_err(x));

//     return 0;
// }

int helper1d::start_test()
{
    int MAX_iteration = 100;
    double epsilon = 1e-8;
    generate_json();
    class Multigrid<1> myMG[4] = {
        Multigrid<1> (MAX_iteration, epsilon, injection, linear, json_path[0]),
        Multigrid<1> (MAX_iteration, epsilon, injection, linear, json_path[1]),
        Multigrid<1> (MAX_iteration, epsilon, injection, linear, json_path[2]),
        Multigrid<1> (MAX_iteration, epsilon, injection, linear, json_path[3])};
    int N_list[4] = {32, 64, 128, 256};
    VectorXd con_V(4), con_FMG(4);
    for(int i = 0 ; i < 4 ; i ++)
    {
        string title = string("$N = ")+std::to_string(N_list[i])+string("$ V-cycle (linear + injection)");
        string filename = string("abserr_N")+std::to_string(N_list[i]);
        VectorXd xV = myMG[i].VC_solver(5,5,true);
        VectorXd xF = myMG[i].FMG_solver(20,20);
        if(i == 3)
        {
            auto vc = myMG[3].get_recorder();
            VectorXd tmp(vc.size());
            for(int i = 0 ; i < vc.size() ; i ++)
                tmp(i) = vc[i];
            output("reduction of relative error for V-cycle ($N=256$)", "VC", "VC", 1, tmp);
        }
        VectorXd eV = get_abs_err(xV);
        VectorXd eF = get_abs_err(xF);
        con_V(i) = eV.lpNorm<Eigen::Infinity>();
        con_FMG(i) = eF.lpNorm<Eigen::Infinity>();
        output(title, filename, "err", 1, eV);
    }
    output("convergence of absolute errors for V-cycle", "rateVC", "convergence", 1, con_V);
    output("convergence of absolute errors for FMG", "rateFMG", "convergence", 1, con_FMG);

    // test epsilon
    auto& mg = myMG[3];
    double eps = 1e-8;
    int1d tmp;
    for(int i = 8 ; i <= 22 ; i +=1)
    {
        mg.set_epsilon(eps);
        mg.VC_solver(5,5,true);
        int t = (mg.get_recorder()).size();
        tmp.push_back(t);
        if(t == 1000)
            break;
        eps /= 10.0;
    }
    VectorXd it_times(tmp.size());
    for(int i = 0 ; i < tmp.size() ; i ++)
        it_times(i) = tmp[i];
    output("the number of iterations", "iter", "iter", 1, it_times);
    // all 8 types
    mg.set_epsilon(1e-8);
    for(int i = 0 ; i < 4 ; i ++)
    {
        mg.set_interpolation(op[i].inter);
        mg.set_restriction(op[i].res);
        VectorXd e = get_abs_err(mg.VC_solver(1,1));
        output("$N=256$ V-cycle "+op[i].name , string("VC")+std::to_string(i), "err", 1, e);
        e = get_abs_err(mg.FMG_solver(20,20));
        output("$N=256$ FMG "+op[i].name , string("FMG")+std::to_string(i), "err", 1, e);
    }
    
    // VectorXd x = myMG[3].VC_solver(1,1,true);
    // double1d v = myMG[3].get_recorder();
    // std::cout << myMG[0].Matrix[0]<< std::endl;
    // std::cout << myMG[0].f<< std::endl;
    // std::cout << myMG[0].LU_solver() << std::endl;
    // std::cout << v.size() << std::endl;
    // std::cout << myMG[0].get_recorder[2] << std::endl;
    // VectorXd e = get_abs_err(x);
    // output("main", "main", "err", 1, e);
    return 0;
}

int helper2d::generate_json()
{
    json_path.clear();
    string prefix = string("../json/") + mark + string("_");
    for(int N = 32 ; N <= 256 ; N *= 2)
    {
        string path = prefix + std::to_string(N) + string(".json");
        json_path.push_back(path);
        class genJson<2> myJson(N, true);
        myJson.set_fval(f);
        if(spectrum != 0)
            myJson.set_sine_domain(spectrum, f, u);
        else if(type == 0) // Dirichlet
            myJson.set_BV0(u);
        else
        {
            myJson.set_BV1(ux, uy, n_vec);
            myJson.set_BV(Vec(0,0),Vec(0,0),u(0,0));
            myJson.set_BV(Vec(0,1),Vec(0,0),u(0,1));
            myJson.set_BV(Vec(1,0),Vec(0,0),u(1,0));
            myJson.set_BV(Vec(1,1),Vec(0,0),u(1,1));
        }
        myJson.write_json(path);
    }
    return 0;
}

VectorXd helper2d::get_abs_err(const VectorXd& _x)
{
    int N = round(std::sqrt(_x.size()));
    assert(N*N == _x.size());
    double h = 1/double(--N);
    VectorXd e(SQUARE(N+1));
    for(int i = 0 ; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
        {
            double x = j*h, y = i*h;
            if(spectrum != 0 && y < spectrum*std::sin(M_PI*x)) // out
                e(to1d(i,j,N)) = 0;
            else
                e(to1d(i,j,N)) = std::abs(u(j*h, i*h) - _x(to1d(i,j,N)));   
        }       
    return e;
}

double helper2d::get_abs_err(const VectorXd& _x, int _norm)
{
    VectorXd e = get_abs_err(_x);
    if(_norm == 0)
        return e.lpNorm<Eigen::Infinity>();
    if(_norm == 1)
        return e.lpNorm<1>();
    if(_norm == 2)
        return e.lpNorm<2>();
    return D_MAX;
}

int helper2d::start_test()
{
    int MAX_iteration = 100;
    double epsilon = 1e-8;
    generate_json();
    // std::cout << 1 << std::endl;
    // int MAX_iteration = 100;
    // double epsilon = 1e-12;
    if(spectrum != 0)
        epsilon = 1e-12;
        
    class Multigrid<2> myMG[4] = {
        Multigrid<2> (MAX_iteration, epsilon, injection, linear, json_path[0]),
        Multigrid<2> (MAX_iteration, epsilon, injection, linear, json_path[1]),
        Multigrid<2> (MAX_iteration, epsilon, injection, linear, json_path[2]),
        Multigrid<2> (MAX_iteration, epsilon, injection, linear, json_path[3])};
    int N_list[4] = {32, 64, 128, 256};
    // std::cout << 2 << std::endl;
    // VectorXd x = myMG[0].VC_solver(1,1);
    // std::cout << get_abs_err(x, 0) << std::endl;
    #define TIME_NOW (std::chrono::high_resolution_clock::now())
    VectorXd conV(4), conF(4);
    double1d tmp_time;
    time_point t1, t2;
    for(int i = 0 ; i <= 3 ; i ++)
    {
        string titleV = string("$N = ")+std::to_string(N_list[i])+string("$ V-cycle (linear + injection)");
        string titleF = string("$N = ")+std::to_string(N_list[i])+string("$ FMG (linear + injection)");
        string filenameV = string("abserr_VC_N")+std::to_string(N_list[i]);
        string filenameF = string("abserr_FMG_N")+std::to_string(N_list[i]);
        
        t1 = TIME_NOW;
        VectorXd xV = myMG[i].VC_solver(5,5,true);
        t2 = TIME_NOW;
        tmp_time.push_back(sec_duration(t1,t2));

        t1 = TIME_NOW;
        VectorXd xF = myMG[i].FMG_solver(20,20);
        t2 = TIME_NOW;
        tmp_time.push_back(sec_duration(t1,t2));

        t1 = TIME_NOW;
        VectorXd xLU = myMG[i].LU_solver();
        t2 = TIME_NOW;
        tmp_time.push_back(sec_duration(t1,t2));

        if(i == 3)
        {
            auto vc = myMG[3].get_recorder();
            VectorXd tmp(vc.size());
            for(int i = 0 ; i < vc.size() ; i ++)
                tmp(i) = vc[i];
            output("reduction of relative error for V-cycle ($N=256$)", "VC", "VC", 2, tmp);
        }
        VectorXd eLU = get_abs_err(xLU);
        VectorXd eV = get_abs_err(xV);
        VectorXd eF = get_abs_err(xF);

        // std::cout << N_list[i] << "_LU: " << eLU.lpNorm<Eigen::Infinity>() << std::endl;
        // std::cout << N_list[i] << "_V: " << eV.lpNorm<Eigen::Infinity>() << std::endl;
        // std::cout << N_list[i] << "_F: " << eF.lpNorm<Eigen::Infinity>() << std::endl;

        conV(i) = eV.lpNorm<Eigen::Infinity>();
        conF(i) = eF.lpNorm<Eigen::Infinity>();
        // // std::cout << "err: " <<  con(i) << std::endl;
        output(titleV, filenameV, "err", 2, eV);
        output(titleF, filenameF, "err", 2, eF);
    }
    // std::cout << "-------------end-------------" << std::endl;
    VectorXd time(tmp_time.size());
    for(int i = 0 ; i < tmp_time.size() ; i ++)
        time(i) = tmp_time[i];
    output("CPU times of different methods", "CPUtime", "time", 2, time);
    output("convergence of absolute errors for V-cycle", "rateVC", "convergence", 2, conV);
    output("convergence of absolute errors for FMG", "rateFMG", "convergence", 2, conF);
    #undef TIME_NOW
    // epsilon
    // test epsilon
    auto& mg = myMG[3];
    double eps = 1e-8;
    int1d tmp;
    for(int i = 8 ; i <= 22 ; i +=1)
    {
        mg.set_epsilon(eps);
        mg.VC_solver(5,5,true);
        int t = (mg.get_recorder()).size();
        tmp.push_back(t);
        if(t >= MAX_iteration)
            break;
        eps /= 10.0;
    }

    VectorXd it_times(tmp.size());
    for(int i = 0 ; i < tmp.size() ; i ++)
        it_times(i) = tmp[i];
    output("the number of iterations", "iter", "iter", 2, it_times);

    //time compare

    
    //time compare
    // auto t1 = std::chrono::steady_clock::now();
    // myMG[2].VC_solver(5,5);
    // auto t2 = std::chrono::steady_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    // std::cout << "time_VC: " << duration.count() << std::endl;
    // t1 = std::chrono::steady_clock::now();
    // myMG[2].LU_solver();
    // t2 = std::chrono::steady_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    // std::cout << "time_LU: " << duration.count() << std::endl;
    return 0;
}