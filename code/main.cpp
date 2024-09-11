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
    int L;
    int n_index = 0;
    int save_more_config = 0;
    std::string folder = "./";
    std::string finfo;
    if (argc == 7)
    {
        std::cout << "Usage: " << argv[0] << " L kappa f g save_more_config folder" << std::endl;
        return 1;

        L = std::atoi(argv[1]);
        Epar.kappa = std::atof(argv[2]);
        Epar.f = std::atof(argv[3]);
        Epar.g = std::atof(argv[4]) * 1.0 / L;
        save_more_config = std::atoi(argv[5]);
        folder = argv[6];

        semiflexible_polymer polymer(L, Epar, 0);
        finfo = "L" + std::string(argv[1]) + "_kappa" + std::string(argv[2]) + "_f" + std::string(argv[3]) + "_gL" + std::string(argv[4]);
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
    }
    else if (argc == 4)
    {
        L = std::atoi(argv[1]);
        n_index = std::atoi(argv[2]);
        folder = argv[3];
        semiflexible_polymer polymer(L, Epar, 1);
        finfo = "L" + std::string(argv[1]) + "_random_run" + std::string(argv[2]);

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
    }
    else
    {
        std::cout << "input error\n";
        return 0;
    }

    std::clock_t c_end = std::clock();
    double time_elapsed = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "Time elapsed: " << time_elapsed << " seconds" << std::endl;

    return 0;
}