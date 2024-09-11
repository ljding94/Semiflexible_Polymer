// Copyright[2024] [Lijie Ding]
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono>
#include <filesystem>
#include "semiflexible_polymer.h"

int main(int argc, char const *argv[])
{
    std::clock_t c_start = std::clock();

    double beta = 1;
    Energy_parameter Epar;

    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm *timeinfo = std::localtime(&now);
    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y%m%d", timeinfo);
    std::string today(buffer);
    std::cout << today << std::endl;

    // precision run with specified parameters
    int L = std::atoi(argv[1]);
    Epar.kappa = std::atof(argv[2]);
    Epar.f = std::atof(argv[3]);
    Epar.g = std::atof(argv[4]) * 1.0 / L;
    int save_more_config = std::atoi(argv[5]);
    std::string folder = argv[6];

    double d_theta = M_PI * 2.0 / 3 / (1.0 + std::sqrt(Epar.kappa));
    std::cout << "d_theta: " << d_theta << std::endl;

    semiflexible_polymer polymer(L, Epar, d_theta);
    std::string finfo = "L" + std::string(argv[1]) + "_kappa" + std::string(argv[2]) + "_f" + std::string(argv[3]) + "_gL" + std::string(argv[4]);

    int bin_num;
    int therm_sweeps;
    int MC_sweeps;
    int step_per_sweep;

    bin_num = 51;
    therm_sweeps = 2000;
    MC_sweeps = 4000;
    step_per_sweep = L * L;

    // polymer.save_polymer_to_file(folder + "/config_" + finfo + "_init.csv");
    // polymer.save_polymer_to_file(folder + "/config_init_" + finfo + ".csv"); // save sample polymer
    polymer.run_simultion(therm_sweeps, MC_sweeps, step_per_sweep, folder, finfo, bin_num, save_more_config);

    std::clock_t c_end = std::clock();
    double time_elapsed = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "Time elapsed: " << time_elapsed << " seconds" << std::endl;

    return 0;
}