//
// Created by useer on 10/17/2018.
//

#ifndef DOUBLE_PENDULUM_TIMING_HPP
#define DOUBLE_PENDULUM_TIMING_HPP

#include <iostream>
#include <chrono>
#include "RKF_algorithm.hpp"
//
// timing codes referenced from Mr. Erik Sevre
//
namespace timing
{
    void timing_RK54_pendulum()
    {
        //
        //set up clock, start timing
        //
        std::chrono::time_point<std::chrono::system_clock> start, stop;
        start = std::chrono::system_clock::now();
        //
        // set up the problem of double pendulum:
        // angles of first and second mass, respectively
        // (0 = vertical downwards, counter-clockwise is positive)
        //
        std::vector<double> theta1;
        std::vector<double> theta2;
        //
        // angle velocities of first and second mass, respectively
        //
        std::vector<double> w1;
        std::vector<double> w2;
        //
        // definition of constant pi
        //
        double pi = acos(-1);
        //
        // Initial conditions. Initial angles of pendulums: theta1=30 degree, theta2=20 degree;
        // Initial angle velocities of pendulums: w1=w2=30 degrees/s
        // unit here is converted into radian
        //
        theta1.insert(theta1.end(), pi*10.0/180.0);
        theta2.insert(theta2.end(), pi*5.0/180.0);
        w1.insert(w1.end(), pi*0.0/180.0);
        w2.insert(w2.end(), pi*0.0/180.0);
        //
        // set up initial time and end time (t0, tf)
        //
        double t0 = 0;
        double tf = 2;
        //
        // call the function to carry out the Runge-Kutta-Fehlberg method
        // (an adaptive Runge-Kutta method of order 5 and 4)
        //
        RK54(t0,tf, theta1, theta2, w1, w2);
        //
        // stop timing
        //
        stop = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = stop - start;
        std::time_t end_time = std::chrono::system_clock::to_time_t(stop);

        std::cout << "finished computation at " << std::ctime(&end_time)
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    }
}
#endif //DOUBLE_PENDULUM_TIMING_HPP
