//
// Created by useer on 10/17/2018.
//

#ifndef DOUBLE_PENDULUM_RKF_ALGORITHM_HPP
#define DOUBLE_PENDULUM_RKF_ALGORITHM_HPP

#include <iostream>
#include <vector>
#include "Runge-Kutta-Fehlberg_functions.hpp"
#include "out_put_data.hpp"
//
// function to check the local error (e) is either equal to or smaller than tolerance (tol)
//
bool err_check(const double &e,
               const double &tol)
{
    return e <= tol;
}
//
// function to find the element that has maximum absolute value in a vector
//
double max(const std::vector<double> &vec)
{
    double result=std::fabs(vec[0]);
    for (auto &v : vec)
    {
        if (result < std::fabs(v)) { result = std::fabs(v); }
    }
    return result;
}
//
// function to carry out the Runge-Kutta-Fehlberg method
// (an adaptive Runge-Kutta method)
//
void RK54(const double &t0,
          const double &tf,
          std::vector<double> theta1,
          std::vector<double> theta2,
          std::vector<double> w1,
          std::vector<double> w2)
{
    //
    // choosing absolute tolerance is 1e-6
    // choosing relative tolerance is 1e-3
    // threshold = absolute tolerance/relative tolerance (after Clive Moler's book)
    //
    double ab_tol = 1e-6;
    double re_tol = 1e-3;
    double threshold = ab_tol / re_tol;
    //
    // create a vector storing time at each step
    //
    std::vector<double> t;
    //
    // insert t0 as the first element of this time vector
    //
    t.push_back(t0);
    //
    // calculate the initial step size h
    // the order "p" of w1, w2, theta1, theta2 (order of y_n+1) at step n+1 is 5, p = 5
    //

    // compute derivatives of angular velocity at t = t0 = 0
    double dw_1=compute_diff_eq_top(theta1[0],theta2[0],w1[0],w2[0]);
    double dw_2=compute_diff_eq_bottom(theta1[0],theta2[0],w1[0],w2[0]);

    // compute initial step size h
    double a = pow(re_tol, 1.0 / 5.0) / max({dw_1/max({w1[0],threshold}),dw_2/max({w2[0],threshold})});
    double h = 0.8 * std::min(std::fabs(t0 - tf), std::fabs(a));

    // set up hmax = 1/10 the time range |tf - t0|
    double hmax = std::fabs(0.1 * (tf - t0));
    //
    // main loop from t0 to tf
    //
    std::cout << "starting main for loop: \n";
    for (size_t i = 0; i < 120000 and t[i] <= tf; ++i) {
        std::cout << "i = " << i << "  t[i] = " << t[i];
        std::cout << "  h = " << h <<std::endl;
        //
        // set up hmin = 16.eps.|t(i)| , with eps = 2^-52
        //
        double hmin = 16 * pow(2, -52) * std::abs(t[i]);
        //
        // adjust so that hmin <= h <= hmax
        //
        if (h > hmax) {
            h = hmax;
        }
        if (h < hmin) {
            (h = hmin);
        }
        //
        // Stretch the step if t is close to tf.
        //
        if (1.1 * h >= (tf - t[i])) {
            h = tf - t[i];
        }

        //
        // compute s1, s2, ..., sk --- see in Runge-Kutta-Fehlberg_functions
        // and return values y_n+1
        //
        theta1.push_back(theta1_new(theta1[i],theta2[i],w1[i],w2[i],h));
        theta2.push_back(theta2_new(theta1[i],theta2[i],w1[i],w2[i],h));
        w1.push_back(w1_new(theta1[i],theta2[i],w1[i],w2[i],h));
        w2.push_back(w2_new(theta1[i],theta2[i],w1[i],w2[i],h));
        //
        // compute local error (le) at t = t_n+1 (i.e. at i+1 step)
        //
        double le_angle = le_theta(theta1[i+1],theta2[i+1],w1[i+1],w2[i+1],h);
        double le_angular_vel = le_w(theta1[i+1],theta2[i+1],w1[i+1],w2[i+1],h);
        //
        // check local error
        //
        // auto check_w = err_check(le_angular_vel,
        //
        // auto check_w = err_check(le_angular_vel,ab_tol);
        //auto b = re_tol * max({w1[i+1],w2[i+1],w1[i],w2[i]});

        //auto c = re_tol * max({theta1[i+1],theta2[i+1],theta1[i],theta2[i]});
        auto le=max({le_angle,le_angular_vel});
        auto max_y_new=max({theta1[i+1],theta2[i+1],w1[i+1],w2[i+1]});
        auto max_y_old=max({theta1[i],theta2[i],w1[i],w2[i]});
        auto check_e = err_check(le, std::max(ab_tol, re_tol * max({max_y_new,max_y_old})));
        //
        // compute new h
        //
        auto h_new = h * 0.8 * pow(ab_tol/le, 1.0 / 5.0);
        //
        // if error <=tolerance (all the check error functions return true)
        // move to the next t until ending the main loop
        //
        //if (check_w and check_theta)
        if (check_e)
        {
            t.push_back(t[i] + h);
            h = h_new;
        } else {
            //
            // if error > tolerance (either one or all of the check error functions return false)
            // redo the main loop with the same t, but the new step size
            //
            w1.pop_back();
            w2.pop_back();
            theta1.pop_back();
            theta2.pop_back();
            h = h_new;
            --i;
        }
        if (std::fabs(t[i]-tf) < 1e-18) {
            break;
        }
    }
    std::cout << "done with main for loop:\n" << std::flush;
    //
    // calculating positions of the center of two masses
    // (x1,y1),(x2,y2) from numerical results
    //
    std::vector<double> x1(t.size());
    std::vector<double> x2(t.size());
    std::vector<double> y1(t.size());
    std::vector<double> y2(t.size());
    pendulum thing;
    for (size_t j = 0; j <t.size() ; ++j) {
        x1[j]=thing.l1*sin(theta1[j]);
        y1[j]=(-1)*thing.l1*cos(theta1[j]);
        x2[j]=x1[j]+ thing.l2*sin(theta2[j]);
        y2[j]=y1[j]- thing.l2*cos(theta2[j]);
    }
    //
    // storing results by calling write vectors function
    //
    std::string filename={"double_pendulum_results.txt"};
    write_vectors_to_file(filename,t,w1,w2,theta1,theta2,x1,y1,x2,y2);

};

//
// end of header file
//
#endif //DOUBLE_PENDULUM_RKF_ALGORITHM_HPP
