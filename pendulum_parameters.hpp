//
// Created by Pham Ngoc Kien on 10/17/2018.
// SNU student ID: 2018-36543
//
//***************************************************************************************************
// This file contains the pendulum parameters (masses and lengths of rods)
//**************************************************************************************************

#ifndef DOUBLE_PENDULUM_PENDULUM_PARAMETERS_HPP
#define DOUBLE_PENDULUM_PENDULUM_PARAMETERS_HPP
struct pendulum {
    // lengths of massless and rigid rods in meters
    double l1 = 1;
    double l2 = 1;

    // masses of two bobs in kg
    double m1 = 1;
    double m2 = 1;
};
#endif //DOUBLE_PENDULUM_PENDULUM_PARAMETERS_HPP
