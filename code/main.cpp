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
    std::string folder;
    double beta = 1;
    Energy_parameter Epar;

    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm *timeinfo = std::localtime(&now);
    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y%m%d", timeinfo);
    std::string today(buffer);
    std::cout << today << std::endl;

    // precision run with specified parameters
    if (argc == 5 || argc == 6)
    {
        int L = std::atoi(argv[1]);
        Epar.kappa = std::atof(argv[2]);
        Epar.f = std::atof(argv[3]);
        Epar.g = std::atof(argv[4]);
        double d_theta = M_PI/(1.0+std::sqrt(Epar.kappa));


        semiflexible_polymer polymer(L, Epar, d_theta);
        std::string finfo = "L" + std::string(argv[1]) + "_kappa" + std::string(argv[2])
                          + "_f" + std::string(argv[3]) + "_g" + std::string(argv[4]);

        int bin_num;
        int therm_sweeps;
        int MC_sweeps;
        int step_per_sweep;
        if(argc == 6)
        {
            // local run
            bin_num = 100;
            therm_sweeps = 1000;
            MC_sweeps = 2000;
            step_per_sweep = L*L*L; //  note: 10L^2 * (1000+1000+2000) for L1000 take 400s
            // use "prog name par* local" for local running
            // used for local running!
            std::cout << "running on local machine\n";
            folder = "../data/scratch_local/" + today;
        } else {
            bin_num = 100;
            therm_sweeps = 2000;
            MC_sweeps = 4000;
            step_per_sweep = 100*L;
            // running on cluster
            std::cout << "running on cluster\n";
            folder = "/home/dingl1/semiflexible_polymer/data_cluster/data_pool"; // dump data to data pool first
        }
        if (!std::filesystem::exists(folder))
            {
                std::cout << today << " folder not exist\n";
                std::cout << "creating folder" << folder << "\n";
                std::filesystem::create_directory(folder);
            }
        // polymer.save_polymer_to_file(folder + "/config_" + finfo + "_init.csv");
        //polymer.save_polymer_to_file(folder + "/config_init_" + finfo + ".csv"); // save sample polymer
        polymer.run_simultion(therm_sweeps, MC_sweeps, step_per_sweep, folder, finfo, bin_num);
    }

    std::clock_t c_end = std::clock();
    double time_elapsed = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "Time elapsed: " << time_elapsed << " seconds" << std::endl;

    return 0;
}