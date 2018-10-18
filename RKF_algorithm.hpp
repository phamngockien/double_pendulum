//
// Created by Pham Ngoc Kien on 10/17/2018.
// SNU student ID: 2018-36543
//
//***********************************************************************************************
// This file produces the Runge-Kutta-Fehlberg algorithm (RK54).
// theta1, theta2, w1, w2 are vectors storing the angles and angular velocities
// of the two pendulum, respectively
// Functions to do the tasks here can be found in Runge-Kutta-Fehlberg_functions header file
//**********************************************************************************************

#ifndef DOUBLE_PENDULUM_RKF_ALGORITHM_HPP
#define DOUBLE_PENDULUM_RKF_ALGORITHM_HPP

#include <iostream>
#include <vector>
#include "Runge-Kutta-Fehlberg_functions.hpp"
#include "out_put_data.hpp"
//*******************************************************************************************
// function to check the local error (e) is either equal to or smaller than tolerance (tol)
//*******************************************************************************************
bool err_check(const double &e,
               const double &tol)
{
    return e <= tol;
}

//**************************************************************************
// function to find the element that has maximum absolute value in a vector
//**************************************************************************
double max(const std::vector<double> &vec)
{
    double result=std::fabs(vec[0]);
    for (auto &v : vec)
    {
        if (result < std::fabs(v)) { result = std::fabs(v); }
    }
    return result;
}

//*************************************************************************************
// function to carry out the Runge-Kutta-Fehlberg method
// (an adaptive Runge-Kutta method of order 5 and 4)
//***********************************************************************************
void RK54(const double &t0,
          const double &tf,
          std::vector<double> theta1,
          std::vector<double> theta2,
          std::vector<double> w1,
          std::vector<double> w2)
{
    //********************************************************************************
    // choosing absolute tolerance, relative tolerance, compute threshold
    // threshold = absolute tolerance/relative tolerance (after Cleve Moler's book)
    //********************************************************************************
    double ab_tol = 1e-6;
    double re_tol = 1e-3;
    double threshold = ab_tol / re_tol;

    // create a vector storing time at each step
    std::vector<double> t;

    // insert t0 as the first element of this time vector
    t.push_back(t0);

    //*********************************************************************************
    // determine initial step size (h)
    //********************************************************************************

    // compute derivatives of angular velocity at t = t0 for the 2 masses
    double dw_1=compute_diff_eq_top(theta1[0],theta2[0],w1[0],w2[0]);
    double dw_2=compute_diff_eq_bottom(theta1[0],theta2[0],w1[0],w2[0]);

    // the order of [w1, w2, theta1, theta2] - (y_n+1) at step n+1 is p = 5
    // compute initial step size h
    double a = pow(re_tol, 1.0 / 5.0) / max({dw_1/max({w1[0],threshold}),dw_2/max({w2[0],threshold})});
    double h = 0.8 * std::min(std::fabs(t0 - tf), std::fabs(a));

    // set hmax = 1/10 the time range |tf - t0|
    double hmax = std::fabs(0.1 * (tf - t0));

    // print out starting main loop
    std::cout << "starting main loop............ \n";
    //******************************************************************
    // main loop from t0 to tf
    //******************************************************************
    for (size_t i = 0; t[i] <= tf; ++i) {
        std::cout << "i = " << i << "  t[i] = " << t[i];
        std::cout << "  h = " << h <<std::endl;

        //*******************************************************************
        // adjust h so that hmin <= h <= hmax
        //*******************************************************************

        // set up hmin = 16.eps.|t(i)| , with eps = 2^-52
        double hmin = 16 * pow(2, -52) * std::abs(t[i]);

        // adjust so that hmin <= h <= hmax
        if (h > hmax) {
            h = hmax;
        }
        if (h < hmin) {
            (h = hmin);
        }

        // Stretch the step if t is close to tf.
        if (1.1 * h >= (tf - t[i])) {
            h = tf - t[i];
        }

        //*******************************************************************************
        // compute s1, s2, ..., sk --- see in Runge-Kutta-Fehlberg_functions
        // and return values y_n+1 for [theta1,theta2,w1,w2] into their container vector
        //********************************************************************************
        theta1.push_back(theta1_new(theta1[i],theta2[i],w1[i],w2[i],h));
        theta2.push_back(theta2_new(theta1[i],theta2[i],w1[i],w2[i],h));
        w1.push_back(w1_new(theta1[i],theta2[i],w1[i],w2[i],h));
        w2.push_back(w2_new(theta1[i],theta2[i],w1[i],w2[i],h));

        //*******************************************************************************
        // compute local error (le) at t = t_n+1 (i.e. at i+1 step)
        // and check the local error
        //*******************************************************************************

        // if you want to calculate both angle and angular velocity error use the underline code
        // double le_angle = le_theta(theta1[i+1],theta2[i+1],w1[i+1],w2[i+1],h);

        // here I compute error with only angular velocities
        double le_angular_vel = le_w(theta1[i+1],theta2[i+1],w1[i+1],w2[i+1],h);

        // if you want to check error with both angle and angular velocity error use the underline code
        //auto le=max({le_angle,le_angular_vel});

         // testing local error with w1,w2 only
         auto le= le_angular_vel;

        // find maximum value of |y_n+1| and |y_n|

        // for both w and theta using these under line codes
        //auto max_y_new=max({theta1[i+1],theta2[i+1],w1[i+1],w2[i+1]});
        //auto max_y_old=max({theta1[i],theta2[i],w1[i],w2[i]});

        // for only using angular velocity (w)
         auto max_y_new=max({w1[i+1],w2[i+1]});
         auto max_y_old=max({w1[i],w2[i]});

        // check local error
        auto check_error = err_check(le, std::max(ab_tol, re_tol * max({max_y_new,max_y_old})));

        //*******************************************************************************
        // compute new h
        //*******************************************************************************
        auto h_new = h * 0.8 * pow(ab_tol/le, 1.0 / 5.0);

        //******************************************************************************
        // the final step in main loop
        //*****************************************************************************
        //
        // if error <=tolerance
        // move to the next time step until ending the main loop
        //
        if (check_error)
        {
            t.push_back(t[i] + h);
            h = h_new;
        } else {
            //
            // if error > tolerance
            // redo the main loop with the same t, but the new step size
            //
            w1.pop_back();
            w2.pop_back();
            theta1.pop_back();
            theta2.pop_back();
            h = h_new;
            --i;
        }

        // break the loop if  | t[i]-tf | < 1e-18
        if (std::fabs(t[i]-tf) < 1e-18) {
            break;
        }
    }
    //*********************************************************************************
    // end of main loop
    //*********************************************************************************
    // print out finish the main loop
    std::cout << "done with main loop \n" << std::flush;

    //*********************************************************************************
    // calculating positions of center of the two masses
    // (x1,y1),(x2,y2) in cartesian coordinate system
    //*********************************************************************************
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
    //*********************************************************************************
    // storing results by calling write vectors function
    //*********************************************************************************
    std::string filename={"double_pendulum_results.txt"};
    write_vectors_to_file(filename,t,w1,w2,theta1,theta2,x1,y1,x2,y2);
};
//
// end of header file
//
#endif //DOUBLE_PENDULUM_RKF_ALGORITHM_HPP
