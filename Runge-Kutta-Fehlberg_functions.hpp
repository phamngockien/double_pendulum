//
// Created by Pham Ngoc Kien on 10/17/2018.
// SNU student ID: 2018-36543
//
//***************************************************************************************************
// This file contains the functions to produce the Runge-Kutta-Fehlberg algorithm (RK54).
// gravity acceleration is 9.8 m/s2
// Pendulum parameters (masses and lengths of rods) can be found in pendulum_parameters header file
//**************************************************************************************************
#ifndef DOUBLE_PENDULUM_RUNGE_KUTTA_FEHLBERG_FUNCTIONS_HPP
#define DOUBLE_PENDULUM_RUNGE_KUTTA_FEHLBERG_FUNCTIONS_HPP

#include "pendulum_parameters.hpp"
#include <cmath>
// ********************************************************************************
// compute differential equations for angular velocities w1 (top) and w2 (bottom)
//********************************************************************************

// compute numerator for the top pendulum to simplify full diff-eq computation
double compute_numerator_top(const double &theta1,
                             const double &theta2,
                             const double &w1,
                             const double &w2,
                             const pendulum &pendulum1=pendulum())
{
    const double g = 9.8;
    return (-1)*g*(2*pendulum1.m1 + pendulum1.m2)*sin(theta1)
           - pendulum1.m2*g*sin(theta1-2*theta2)
           - 2*sin(theta1-theta2)*pendulum1.m2*(w2*w2*pendulum1.l2+w1*w1*pendulum1.l1*cos(theta1-theta2));
}

// compute denominator for the top pendulum to simplify full diff-eq computation
double compute_denominator_top(const double &theta1,
                               const double &theta2,
                               const pendulum &pendulum1=pendulum())
{
    return pendulum1.l1*(2*pendulum1.m1 + pendulum1.m2 - pendulum1.m2*cos(2*theta1-2*theta2));
}

//  compute solution to the diff eq for the top pendulum
double compute_diff_eq_top(const double &theta1,
                           const double &theta2,
                           const double &w1,
                           const double &w2)
{
    return compute_numerator_top(theta1,theta2,w1,w2)/compute_denominator_top(theta1,theta2);
}
//---------------------------------------------------------------------------------------------------------------
// compute numerator for the bottom pendulum to simplify full diff-eq computation
double compute_numerator_bottom(const double &theta1,
                                const double &theta2,
                                const double &w1,
                                const double &w2,
                                const pendulum &pendulum1=pendulum())
{
    const double g = 9.8;
    return 2*sin(theta1-theta2)
           * (w1*w1*pendulum1.l1*(pendulum1.m1 + pendulum1.m2)+g*(pendulum1.m1 + pendulum1.m2)*cos(theta1)
              + w2*w2*pendulum1.l2*pendulum1.m2*cos(theta1-theta1) );
}

// compute denominator for the bottom pendulum to simplify full diff-eq computation
double compute_denominator_bottom(const double &theta1,
                                  const double &theta2,
                                  const pendulum &pendulum1=pendulum())
{
    return pendulum1.l2*(2*pendulum1.m1 + pendulum1.m2 - pendulum1.m2*cos(2*theta1-2*theta2));
}

//  compute solution to the diff eq for the bottom pendulum
double compute_diff_eq_bottom(const double &theta1,
                              const double &theta2,
                              const double &w1,
                              const double &w2)
{
    return compute_numerator_bottom(theta1,theta2,w1,w2,pendulum())/compute_denominator_bottom(theta1,theta2,pendulum());
}

//********************************************************************************************
// compute s1, s2, ..., s6 in the Runge-Kutta-Fehlberg method(RK54).
//********************************************************************************************

//-------------------------------------------------------------------
// calculate s1 for theta1, theta2, w1, w2, respectively
//-------------------------------------------------------------------
double s1_theta1(const double &w1)
{
    return w1;
}
//-------------------------------------------------------------------
double s1_theta2(const double &w2)
{
    return w2;
}
//-------------------------------------------------------------------
double s1_w1(const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    return compute_diff_eq_top(theta1,theta2,w1,w2);
}
//-------------------------------------------------------------------
double s1_w2(const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    return compute_diff_eq_bottom(theta1,theta2,w1,w2);
}

//-------------------------------------------------------------------
// calculate s2 for theta1, theta2, w1, w2, respectively
//-------------------------------------------------------------------
double s2_theta1(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h)
{
    return s1_theta1(w1) + s1_w1(theta1,theta2,w1,w2)*(1.0/4.0)*h;
}
//-------------------------------------------------------------------
double s2_theta2(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h)
{
    return s1_theta2(w2) + s1_w2(theta1,theta2,w1,w2)*(1.0/4.0)*h;
}
//-------------------------------------------------------------------
double s2_w1(const double &h,
             const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    double theta1_2 = theta1 + s1_theta1(w1)*(1.0/4.0)*h;
    double theta2_2 = theta2 + s1_theta2(w2)*(1.0/4.0)*h;
    double w1_2= w1 + s1_w1(theta1,theta2,w1,w2)*(1.0/4.0)*h;
    double w2_2= w2 + s1_w2(theta1,theta2,w1,w2)*(1.0/4.0)*h;
    return compute_diff_eq_top(theta1_2,theta2_2,w1_2,w2_2);
}
//-------------------------------------------------------------------
double s2_w2(const double &h,
             const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    double theta1_2 = theta1 + s1_theta1(w1)*(1.0/4.0)*h;
    double theta2_2 = theta2 + s1_theta2(w2)*(1.0/4.0)*h;
    double w1_2= w1 + s1_w1(theta1,theta2,w1,w2)*(1.0/4.0)*h;
    double w2_2= w2 + s1_w2(theta1,theta2,w1,w2)*(1.0/4.0)*h;
    return compute_diff_eq_bottom(theta1_2,theta2_2,w1_2,w2_2);
}

//-------------------------------------------------------------------
// calculate s3 for theta1, theta2, w1, w2, respectively
//-------------------------------------------------------------------
double s3_theta1(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h)
{
    return s1_theta1(w1) + s1_w1(theta1,theta2,w1,w2)*(3.0/32.0)*h
           + s2_w1(h,theta1,theta2,w1,w2)*(9.0/32.0)*h;
}
//-------------------------------------------------------------------
double s3_theta2(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h)
{
    return s1_theta2(w2) + s1_w2(theta1,theta2,w1,w2)*(3.0/32.0)*h
           + s2_w2(h,theta1,theta2,w1,w2)*(9.0/32.0)*h;
}
//-------------------------------------------------------------------
double s3_w1(const double &h,
             const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    double theta1_2 = theta1 + s1_theta1(w1)*(3.0/32.0)*h
                      + s2_theta1(theta1,theta2,w1,w2,h)*(9.0/32.0)*h;
    double theta2_2 = theta2 + s1_theta2(w2)*(3.0/32.0)*h
                      + s2_theta2(theta1,theta2,w1,w2,h)*(9.0/32.0)*h;
    double w1_2= w1 + s1_w1(theta1,theta2,w1,w2)*(3.0/32.0)*h
                 + s2_w1(h,theta1,theta2,w1,w2)*(9.0/32.0)*h ;
    double w2_2= w2 + s1_w2(theta1,theta2,w1,w2)*(3.0/32.0)*h
                 + s2_w2(h,theta1,theta2,w1,w2)*(9.0/32.0)*h;
    return compute_diff_eq_top(theta1_2,theta2_2,w1_2,w2_2);
}
//-------------------------------------------------------------------
double s3_w2(const double &h,
             const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    double theta1_2 = theta1 + s1_theta1(w1)*(3.0/32.0)*h
                      + s2_theta1(theta1,theta2,w1,w2,h)*(9.0/32.0)*h;
    double theta2_2 = theta2 + s1_theta2(w2)*(3.0/32.0)*h
                      + s2_theta2(theta1,theta2,w1,w2,h)*(9.0/32.0)*h;
    double w1_2= w1 + s1_w1(theta1,theta2,w1,w2)*(3.0/32.0)*h
                 + s2_w1(h,theta1,theta2,w1,w2)*(9.0/32.0)*h ;
    double w2_2= w2 + s1_w2(theta1,theta2,w1,w2)*(3.0/32.0)*h
                 + s2_w2(h,theta1,theta2,w1,w2)*(9.0/32.0)*h;
    return compute_diff_eq_bottom(theta1_2,theta2_2,w1_2,w2_2);
}

//-------------------------------------------------------------------
// calculate s4 for theta1, theta2, w1, w2, respectively
//-------------------------------------------------------------------
double s4_theta1(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h)
{
    return s1_theta1(w1) + s1_w1(theta1,theta2,w1,w2)*(1932.0/2197.0)*h
           + s2_w1(h,theta1,theta2,w1,w2)*(-7200.0/2197.0)*h
           + s3_w1(h,theta1,theta2,w1,w2)*(7296.0/2197.0)*h;
}
//-------------------------------------------------------------------
double s4_theta2(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h)
{
    return s1_theta2(w2) + s1_w2(theta1,theta2,w1,w2)*(1932.0/2197.0)*h
           + s2_w2(h,theta1,theta2,w1,w2)*(-7200.0/2197.0)*h
           + s3_w2(h,theta1,theta2,w1,w2)*(7296.0/2197.0)*h;
}
//-------------------------------------------------------------------
double s4_w1(const double &h,
             const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    double theta1_2 = theta1 + s1_theta1(w1)*(1932.0/2197.0)*h
                      + s2_theta1(theta1,theta2,w1,w2,h)*(-7200.0/2197.0)*h
                      + s3_theta1(theta1,theta2,w1,w2,h)*(7296.0/2197.0)*h;
    double theta2_2 = theta2 + s1_theta2(w2)*(1932.0/2197.0)*h
                      + s2_theta2(theta1,theta2,w1,w2,h)*(-7200.0/2197.0)*h
                      + s3_theta2(theta1,theta2,w1,w2,h)*(7296.0/2197.0)*h;
    double w1_2= w1 + s1_w1(theta1,theta2,w1,w2)*(1932.0/2197.0)*h
                 + s2_w1(h,theta1,theta2,w1,w2)*(-7200.0/2197.0)*h
                 + s3_w1(h,theta1,theta2,w1,w2)*(7296.0/2197.0)*h ;
    double w2_2= w2 + s1_w2(theta1,theta2,w1,w2)*(1932.0/2197.0)*h
                 + s2_w2(h,theta1,theta2,w1,w2)*(-7200.0/2197.0)*h
                 + s3_w2(h,theta1,theta2,w1,w2)*(7296.0/2197.0)*h ;
    return compute_diff_eq_top(theta1_2,theta2_2,w1_2,w2_2);
}
//-------------------------------------------------------------------
double s4_w2(const double &h,
             const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    double theta1_2 = theta1 + s1_theta1(w1)*(1932.0/2197.0)*h
                      + s2_theta1(theta1,theta2,w1,w2,h)*(-7200.0/2197.0)*h
                      + s3_theta1(theta1,theta2,w1,w2,h)*(7296.0/2197.0)*h;
    double theta2_2 = theta2 + s1_theta2(w2)*(1932.0/2197.0)*h
                      + s2_theta2(theta1,theta2,w1,w2,h)*(-7200.0/2197.0)*h
                      + s3_theta2(theta1,theta2,w1,w2,h)*(7296.0/2197.0)*h;
    double w1_2= w1 + s1_w1(theta1,theta2,w1,w2)*(1932.0/2197.0)*h
                 + s2_w1(h,theta1,theta2,w1,w2)*(-7200.0/2197.0)*h
                 + s3_w1(h,theta1,theta2,w1,w2)*(7296.0/2197.0)*h ;
    double w2_2= w2 + s1_w2(theta1,theta2,w1,w2)*(1932.0/2197.0)*h
                 + s2_w2(h,theta1,theta2,w1,w2)*(-7200.0/2197.0)*h
                 + s3_w2(h,theta1,theta2,w1,w2)*(7296.0/2197.0)*h ;
    return compute_diff_eq_bottom(theta1_2,theta2_2,w1_2,w2_2);
}
//-------------------------------------------------------------------
// calculate s5 for theta1, theta2, w1, w2, respectively
//-------------------------------------------------------------------
double s5_theta1(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h)
{
    return s1_theta1(w1) + s1_w1(theta1,theta2,w1,w2)*(439.0/216.0)*h
           + s2_w1(h,theta1,theta2,w1,w2)*(-8.0)*h
           + s3_w1(h,theta1,theta2,w1,w2)*(3680.0/513.0)*h
           + s4_w1(h,theta1,theta2,w1,w2)*(-845.0/4104)*h;
}
//-------------------------------------------------------------------
double s5_theta2(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h)
{
    return s1_theta2(w2) + s1_w2(theta1,theta2,w1,w2)*(439.0/216.0)*h
           + s2_w2(h,theta1,theta2,w1,w2)*(-8.0)*h
           + s3_w2(h,theta1,theta2,w1,w2)*(3680.0/513.0)*h
           + s4_w2(h,theta1,theta2,w1,w2)*(-845.0/4104)*h;
}
//-------------------------------------------------------------------
double s5_w1(const double &h,
             const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    double theta1_2 = theta1 + s1_theta1(w1)*(439.0/216.0)*h
                      + s2_theta1(theta1,theta2,w1,w2,h)*(-8.0)*h
                      + s3_theta1(theta1,theta2,w1,w2,h)*(3680.0/513.0)*h
                      + s4_theta1(theta1,theta2,w1,w2,h)*(-845.0/4104)*h;
    double theta2_2 = theta2 + s1_theta2(w2)*(439.0/216.0)*h
                      + s2_theta2(theta1,theta2,w1,w2,h)*(-8.0)*h
                      + s3_theta2(theta1,theta2,w1,w2,h)*(3680.0/513.0)*h
                      + s4_theta2(theta1,theta2,w1,w2,h)*(-845.0/4104)*h;
    double w1_2= w1 + s1_w1(theta1,theta2,w1,w2)*(439.0/216.0)*h
                 + s2_w1(h,theta1,theta2,w1,w2)*(-8.0)*h
                 + s3_w1(h,theta1,theta2,w1,w2)*(3680.0/513.0)*h
                 + s4_w1(h,theta1,theta2,w1,w2)*(-845.0/4104)*h;
    double w2_2= w2 + s1_w2(theta1,theta2,w1,w2)*(439.0/216.0)*h
                 + s2_w2(h,theta1,theta2,w1,w2)*(-8.0)*h
                 + s3_w2(h,theta1,theta2,w1,w2)*(3680.0/513.0)*h
                 + s4_w2(h,theta1,theta2,w1,w2)*(-845.0/4104)*h;
    return compute_diff_eq_top(theta1_2,theta2_2,w1_2,w2_2);
}
//-------------------------------------------------------------------
double s5_w2(const double &h,
             const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    double theta1_2 = theta1 + s1_theta1(w1)*(439.0/216.0)*h
                      + s2_theta1(theta1,theta2,w1,w2,h)*(-8.0)*h
                      + s3_theta1(theta1,theta2,w1,w2,h)*(3680.0/513.0)*h
                      + s4_theta1(theta1,theta2,w1,w2,h)*(-845.0/4104)*h;
    double theta2_2 = theta2 + s1_theta2(w2)*(439.0/216.0)*h
                      + s2_theta2(theta1,theta2,w1,w2,h)*(-8.0)*h
                      + s3_theta2(theta1,theta2,w1,w2,h)*(3680.0/513.0)*h
                      + s4_theta2(theta1,theta2,w1,w2,h)*(-845.0/4104)*h;
    double w1_2= w1 + s1_w1(theta1,theta2,w1,w2)*(439.0/216.0)*h
                 + s2_w1(h,theta1,theta2,w1,w2)*(-8.0)*h
                 + s3_w1(h,theta1,theta2,w1,w2)*(3680.0/513.0)*h
                 + s4_w1(h,theta1,theta2,w1,w2)*(-845.0/4104)*h;
    double w2_2= w2 + s1_w2(theta1,theta2,w1,w2)*(439.0/216.0)*h
                 + s2_w2(h,theta1,theta2,w1,w2)*(-8.0)*h
                 + s3_w2(h,theta1,theta2,w1,w2)*(3680.0/513.0)*h
                 + s4_w2(h,theta1,theta2,w1,w2)*(-845.0/4104)*h;
    return compute_diff_eq_bottom(theta1_2,theta2_2,w1_2,w2_2);
}
//-------------------------------------------------------------------
// calculate s6 for theta1, theta2, w1, w2, respectively
//-------------------------------------------------------------------
double s6_theta1(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h)
{
    return s1_theta1(w1) + s1_w1(theta1,theta2,w1,w2)*(-8.0/27.0)*h
           + s2_w1(h,theta1,theta2,w1,w2)*(2.0)*h
           + s3_w1(h,theta1,theta2,w1,w2)*(-3544.0/2565.0)*h
           + s4_w1(h,theta1,theta2,w1,w2)*(1859.0/4104.0)*h
           + s5_w1(h,theta1,theta2,w1,w2)*(-11.0/40.0)*h;
}
//-------------------------------------------------------------------
double s6_theta2(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h)
{
    return s1_theta2(w2) + s1_w2(theta1,theta2,w1,w2)*(-8.0/27.0)*h
           + s2_w2(h,theta1,theta2,w1,w2)*(2.0)*h
           + s3_w2(h,theta1,theta2,w1,w2)*(-3544.0/2565.0)*h
           + s4_w2(h,theta1,theta2,w1,w2)*(1859.0/4104.0)*h
           + s5_w2(h,theta1,theta2,w1,w2)*(-11.0/40.0)*h;
}
//-------------------------------------------------------------------
double s6_w1(const double &h,
             const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    double theta1_2 = theta1 + s1_theta1(w1)*(-8.0/27.0)*h
                      + s2_theta1(theta1,theta2,w1,w2,h)*(2.0)*h
                      + s3_theta1(theta1,theta2,w1,w2,h)*(-3544.0/2565.0)*h
                      + s4_theta1(theta1,theta2,w1,w2,h)*(1859.0/4104.0)*h
                      + s5_theta1(theta1,theta2,w1,w2,h)*(-11.0/40.0)*h;
    double theta2_2 = theta2 + s1_theta2(w1)*(-8.0/27.0)*h
                      + s2_theta2(theta1,theta2,w1,w2,h)*(2.0)*h
                      + s3_theta2(theta1,theta2,w1,w2,h)*(-3544.0/2565.0)*h
                      + s4_theta2(theta1,theta2,w1,w2,h)*(1859.0/4104.0)*h
                      + s5_theta2(theta1,theta2,w1,w2,h)*(-11.0/40.0)*h;
    double w1_2= w1 + s1_w1(theta1,theta2,w1,w2)*(-8.0/27.0)*h
                 + s2_w1(h,theta1,theta2,w1,w2)*(2.0)*h
                 + s3_w1(h,theta1,theta2,w1,w2)*(-3544.0/2565.0)*h
                 + s4_w1(h,theta1,theta2,w1,w2)*(1859.0/4104.0)*h
                 + s5_w1(h,theta1,theta2,w1,w2)*(-11.0/40.0)*h;
    double w2_2= w2 + s1_w2(theta1,theta2,w1,w2)*(-8.0/27.0)*h
                 + s2_w2(h,theta1,theta2,w1,w2)*(2.0)*h
                 + s3_w2(h,theta1,theta2,w1,w2)*(-3544.0/2565.0)*h
                 + s4_w2(h,theta1,theta2,w1,w2)*(1859.0/4104.0)*h
                 + s5_w2(h,theta1,theta2,w1,w2)*(-11.0/40.0)*h;
    return compute_diff_eq_top(theta1_2,theta2_2,w1_2,w2_2);
}
//-------------------------------------------------------------------
double s6_w2(const double &h,
             const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2)
{
    double theta1_2 = theta1 + s1_theta1(w1)*(-8.0/27.0)*h
                      + s2_theta1(theta1,theta2,w1,w2,h)*(2.0)*h
                      + s3_theta1(theta1,theta2,w1,w2,h)*(-3544.0/2565.0)*h
                      + s4_theta1(theta1,theta2,w1,w2,h)*(1859.0/4104.0)*h
                      + s5_theta1(theta1,theta2,w1,w2,h)*(-11.0/40.0)*h;
    double theta2_2 = theta2 + s1_theta2(w1)*(-8.0/27.0)*h
                      + s2_theta2(theta1,theta2,w1,w2,h)*(2.0)*h
                      + s3_theta2(theta1,theta2,w1,w2,h)*(-3544.0/2565.0)*h
                      + s4_theta2(theta1,theta2,w1,w2,h)*(1859.0/4104.0)*h
                      + s5_theta2(theta1,theta2,w1,w2,h)*(-11.0/40.0)*h;
    double w1_2= w1 + s1_w1(theta1,theta2,w1,w2)*(-8.0/27.0)*h
                 + s2_w1(h,theta1,theta2,w1,w2)*(2.0)*h
                 + s3_w1(h,theta1,theta2,w1,w2)*(-3544.0/2565.0)*h
                 + s4_w1(h,theta1,theta2,w1,w2)*(1859.0/4104.0)*h
                 + s5_w1(h,theta1,theta2,w1,w2)*(-11.0/40.0)*h;
    double w2_2= w2 + s1_w2(theta1,theta2,w1,w2)*(-8.0/27.0)*h
                 + s2_w2(h,theta1,theta2,w1,w2)*(2.0)*h
                 + s3_w2(h,theta1,theta2,w1,w2)*(-3544.0/2565.0)*h
                 + s4_w2(h,theta1,theta2,w1,w2)*(1859.0/4104.0)*h
                 + s5_w2(h,theta1,theta2,w1,w2)*(-11.0/40.0)*h;
    return compute_diff_eq_bottom(theta1_2,theta2_2,w1_2,w2_2);
}
//**************************************************************************************
// calculate y (t = t_n+1) for theta1, theta2, w1, w2
//**************************************************************************************
//
double theta1_new(const double &theta1,
                  const double &theta2,
                  const double &w1,
                  const double &w2,
                  const double &h)
{
    return theta1 + h * ( (16.0/135.0)*s1_theta1(w1) + (0.0)*s2_theta1(theta1,theta2,w1,w2,h)
                          + (6656.0/12825.0)*s3_theta1(theta1,theta2,w1,w2,h)
                          + (28561.0/56430.0)*s4_theta1(theta1,theta2,w1,w2,h)
                          + (-9.0/50.0)*s5_theta1(theta1,theta2,w1,w2,h)
                          + (2.0/55.0)*s6_theta1(theta1,theta2,w1,w2,h) );
}
//-------------------------------------------------------------------
double theta2_new(const double &theta1,
                  const double &theta2,
                  const double &w1,
                  const double &w2,
                  const double &h)
{
    return theta2 + h * ( (16.0/135.0)*s1_theta2(w2) + (0.0)*s2_theta2(theta1,theta2,w1,w2,h)
                          + (6656.0/12825.0)*s3_theta2(theta1,theta2,w1,w2,h)
                          + (28561.0/56430.0)*s4_theta2(theta1,theta2,w1,w2,h)
                          + (-9.0/50.0)*s5_theta2(theta1,theta2,w1,w2,h)
                          + (2.0/55.0)*s6_theta2(theta1,theta2,w1,w2,h) );
}
//-------------------------------------------------------------------
double w1_new(const double &theta1,
              const double &theta2,
              const double &w1,
              const double &w2,
              const double &h)
{
    return w1 + h * ( (16.0/135.0)*s1_w1(theta1,theta2,w1,w2)
                      + (0.0)*s2_w1(h,theta1,theta2,w1,w2)
                      + (6656.0/12825.0)*s3_w1(h,theta1,theta2,w1,w2)
                      + (28561.0/56430.0)*s4_w1(h,theta1,theta2,w1,w2)
                      + (-9.0/50.0)*s5_w1(h,theta1,theta2,w1,w2)
                      + (2.0/55.0)*s6_w1(h,theta1,theta2,w1,w2) );
}
//-------------------------------------------------------------------
double w2_new(const double &theta1,
              const double &theta2,
              const double &w1,
              const double &w2,
              const double &h)
{
    return w2 + h * ( (16.0/135.0)*s1_w2(theta1,theta2,w1,w2)
                      + (0.0)*s2_w2(h,theta1,theta2,w1,w2)
                      + (6656.0/12825.0)*s3_w2(h,theta1,theta2,w1,w2)
                      + (28561.0/56430.0)*s4_w2(h,theta1,theta2,w1,w2)
                      + (-9.0/50.0)*s5_w2(h,theta1,theta2,w1,w2)
                      + (2.0/55.0)*s6_w2(h,theta1,theta2,w1,w2) );
}
//**************************************************************************************
// compute z_n+1 for theta1, theta2, w1, w2
//**************************************************************************************
double theta1_z_new(const double &theta1,
                  const double &theta2,
                  const double &w1,
                  const double &w2,
                  const double &h)
{
    return theta1 + h * ( (25.0/216.0) * s1_theta1(w1)
                               + (0.0) * s2_theta1(theta1, theta2, w1, w2, h)
                     + (1408.0/2565.0) * s3_theta1(theta1, theta2, w1, w2, h)
                     + (2197.0/4104.0) * s4_theta1(theta1, theta2, w1, w2, h)
                         + (- 1.0/5.0) * s5_theta1(theta1, theta2, w1, w2, h) );
}
//-------------------------------------------------------------------
double theta2_z_new(const double &theta1,
                  const double &theta2,
                  const double &w1,
                  const double &w2,
                  const double &h)
{
    return theta2 + h * ( (25.0/216.0) * s1_theta2(w2)
                               + (0.0) * s2_theta2(theta1, theta2, w1, w2, h)
                     + (1408.0/2565.0) * s3_theta2(theta1, theta2, w1, w2, h)
                     + (2197.0/4104.0) * s4_theta2(theta1, theta2, w1, w2, h)
                         + (- 1.0/5.0) * s5_theta2(theta1, theta2, w1, w2, h) );
}
//-------------------------------------------------------------------
double w1_z_new(const double &theta1,
              const double &theta2,
              const double &w1,
              const double &w2,
              const double &h)
{
    return w1 +h * ( (25.0/216.0) * s1_w1(theta1,theta2,w1,w2)
                          + (0.0) * s2_w1(h,theta1,theta2,w1,w2)
                + (1408.0/2565.0) * s3_w1(h,theta1,theta2,w1,w2)
                + (2197.0/4104.0) * s4_w1(h,theta1,theta2,w1,w2)
                    + (- 1.0/5.0) * s5_w1(h,theta1,theta2,w1,w2) );
}
//-------------------------------------------------------------------
double w2_z_new(const double &theta1,
              const double &theta2,
              const double &w1,
              const double &w2,
              const double &h)
{
    return w2 + h * ( (25.0/216.0) * s1_w2(theta1,theta2,w1,w2)
                           + (0.0) * s2_w2(h,theta1,theta2,w1,w2)
                 + (1408.0/2565.0) * s3_w2(h,theta1,theta2,w1,w2)
                 + (2197.0/4104.0) * s4_w2(h,theta1,theta2,w1,w2)
                     + (- 1.0/5.0) * s5_w2(h,theta1,theta2,w1,w2) );
}
//****************************************************************************
// calculate the error e_n+1 = y_n+1 - z_n+1
//****************************************************************************
//------------------------------------------------------
// local error for angles (le_theta1, le_theta2)
//------------------------------------------------------
double le_theta1(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h) {
    return theta1_new(theta1,theta2,w1,w2,h)  - theta1_z_new(theta1,theta2,w1,w2,h);
}
//-----------------------------------------------------------------
double le_theta2(const double &theta1,
                 const double &theta2,
                 const double &w1,
                 const double &w2,
                 const double &h) {
    return theta2_new(theta1,theta2,w1,w2,h)  - theta2_z_new(theta1,theta2,w1,w2,h);
}
//-----------------------------------------------------------------
// maximum absolute value of local errors of angles
double le_theta (const double &le_theta1,
                 const double &le_theta2) {
    return std::fmax(std::fabs(le_theta1),std::fabs(le_theta2));
}
//------------------------------------------------------
// local error for angular velocities (le_w1, le_w2)
//------------------------------------------------------
double le_w1 (const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2,
             const double &h){
    return w1_new(theta1,theta2,w1,w2,h) - w1_z_new(theta1,theta2,w1,w2,h);
}
//-----------------------------------------------------------------
double le_w2 (const double &theta1,
             const double &theta2,
             const double &w1,
             const double &w2,
             const double &h){
    return w2_new(theta1,theta2,w1,w2,h)-  w2_z_new(theta1,theta2,w1,w2,h);
}
//-----------------------------------------------------------------
// maximum absolute value of local errors of angular velocities
double le_w (const double &le_w1,
             const double &le_w2)
{
    return std::fmax(std::fabs(le_w1),std::fabs(le_w2));
}
#endif //DOUBLE_PENDULUM_RUNGE_KUTTA_FEHLBERG_FUNCTIONS_HPP
