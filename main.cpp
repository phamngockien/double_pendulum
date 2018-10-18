//
// Created by Pham Ngoc Kien on 10/17/2018.
// SNU student ID: 2018-36543
//
//***************************************************************************************************
// This file contains a main function to produce the ODE solver by Runge-Kutta-Fehlberg method
// Algorithm, calculating and timing functions can be seen in supporting header files.
// Pendulum parameters (masses and lengths of rods) can be found in pendulum_parameters header file
//**************************************************************************************************
#include <iostream>
#include "timing.hpp"

int main() {
    timing::timing_RK54_pendulum();
}