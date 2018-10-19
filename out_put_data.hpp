//
// Created by useer on 10/17/2018.
//

#ifndef DOUBLE_PENDULUM_OUT_PUT_DATA_HPP
#define DOUBLE_PENDULUM_OUT_PUT_DATA_HPP
#include <string>
#include <fstream>
#include <vector>

//
// function to write vector to file
// columns of this file containing:
// time (s),
// angular velocites w1, w2 (rad/s),
// angles theta1,theta2 (rad),
// pendulum positions x1,y1,x2,y2 (m),
// errors of angles and angular velocity
void write_vectors_to_file(
        const std::string &filename,
        const std::vector<double> &t,
        const std::vector<double> &w1,
        const std::vector<double> &w2,
        const std::vector<double> &theta1,
        const std::vector<double> &theta2,
        const std::vector<double> &x1,
        const std::vector<double> &y1,
        const std::vector<double> &x2,
        const std::vector<double> &y2,
        const std::vector<double> &w1_error,
        const std::vector<double> &w2_error,
        const std::vector<double> &theta1_error,
        const std::vector<double> &theta2_error)
{
    std::ofstream outfile(filename);
    for (size_t i = 0; i < t.size(); ++i) {
        outfile << t[i] << "\t" << w1[i] << "\t" << w2[i] << "\t"
                << theta1[i] << "\t" << theta2[i] << "\t"
                << x1[i] << "\t" << y1[i] << "\t"
                << x2[i] << "\t" << y2[i] << "\t"
                << w1_error[i] << "\t" << w2_error[i] << "\t"
                << theta1_error[i] << "\t" << theta2_error[i] << "\n";
    }
    outfile.close();
}
#endif //DOUBLE_PENDULUM_OUT_PUT_DATA_HPP
