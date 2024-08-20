#ifndef _SEMIFLEXIBLE_POLYMER_H
#define _SEMIFLEXIBLE_POLYMER_H
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

// hamiltonion parameters
struct Energy_parameter
{
    double kappa = 0; // bending modulii
    double f = 0;     // force term, here assume in x direction
    double g = 0;     // shear term, here assume in z direction

    // Energy:
    // E = sum {\kappa/2*(t_i - t_{i-1})^2 - (f*t_i(x) + g* |r_i(z)|*t_i(x)}
    // abs z |r_i(z)| due to shear streching direction alon x switch sign at z=0
};
struct observable
{
    // put all observable interested here
    double E;                            // energy
    double Tb;                           // total bending: bending energy = 0.5*kappa*Tb
    double X;                            // end to end X distance  |(r(L-1) + t(L-1) - r(0)).x|
    double Y;                            // end to end Y distance  |(r(L-1) + t(L-1) - r(0)).y|
    double Z;                            // end to end Z distance  |(r(L-1) + t(L-1) - r(0)).z|
    double XsignZ;                       // end to end X distance  times sign of Z, to eliminate +- symmetry
    double ZsignX;                       // end to end Z distance  times sign of X, to eliminate +- symmetry
    double R;                            // end to end distance  |r(L-1) + t(L-1) - r(0)|
    double R2;                           // end to end distance square
    double Rg2;                           // radius of gyration
    double Sxx, Syy, Szz, Sxy, Sxz, Syz; // gyration tensor components

    std::vector<double> Sq{};  // structure factor
    std::vector<double> qB{};  // qB for Sq, B is the bead-bead distance, which is 1 in our unit
    std::vector<double> tts{}; // tangent-tangent correlation function, versus contour distance
    std::vector<double> spB{}; // s/B for tts calculation
};
struct bead
{
    // configuration related
    std::vector<double> r{0, 0, 0};     // position (x,y,z)
    std::vector<double> t{1.0, 0, 0.0}; // tangent, point to next bead R' = R+t, bead-bead distance B is unit vector
};

class semiflexible_polymer
{
public:
    double beta;    // system temperature
    int L;          // length, (L+1) beads in total  (0)--(1)--(2)-- ... --(L-2)--(L-1)--(L)
    int fix_bead_0; // fix the 0th bead

    Energy_parameter Epar;
    std::vector<bead> polymer; // the polymer
    double E_sys;              // system energy
    double Tb_sys;             // system bending

    double d_theta; // for rot update

    // randomnumber generators
    std::mt19937 gen;
    std::uniform_real_distribution<> rand_uni; // uniform distribution
    std::normal_distribution<> rand_norm;      // normal distribution

    std::vector<double> rand_uni_vec(); // generate a random unit vector

    // initialization
    semiflexible_polymer(double L_, Energy_parameter Epar_, double d_theta_, int fixe_bead_0_ = 1);
    // MCMC update
    int update_bead_concerted_rotation(int bead_i, int bead_j);
    // con-rot update in a corn with angle d_theta
    int update_bead_tangent_rotation(int bead_i);
    // rotage the tangent of bead i, every bead > i position will be changed

    // check self-avoid condition
    int satisfy_self_avoiding_condition(int bead_i);                      // check if polymer[i] satisfy self-avoiding condition with all other beads
    int satisfy_self_avoiding_condition_by_group(int bead_i, int bead_j); //
    // check if polymer[:i] and polymer[i:] satisfy self-avoiding condition with all other beads

    // some observable measurement
    observable measure_observable(int bin_num);
    // pair distribution function
    std::vector<double> calc_structure_factor(std::vector<double> qB);          // structure factor of the polymer, orientational averaged: ref: Pedersen 1996 equ 1
    std::vector<double> calc_rod_structure_factor(std::vector<double> qB);      // calculate the structure factor of a rod
    std::vector<double> calc_tangent_pair_correlation(std::vector<double> spB); // calculate the pair tangent correlation distribution function
    std::vector<double> calc_gyration_tensor();                                 // calculate the gyration tensor
    double calc_radius_of_gyration_square();                                           // radius of gyration Rg^2

    // experiment
    void save_polymer_to_file(std::string filename);
    void save_observable_to_file(std::string filename, std::vector<observable> obs_ensemble, bool save_detail = false);

    void run_simultion(int therm_sweep, int MC_sweeps, int step_per_sweep, std::string folder, std::string finfo, int bin_num, int save_more_config);

    // spB stand for s per B or s/P

    // little tool
    std::vector<double> cross_product(std::vector<double> a, std::vector<double> b); // cross product of a and b
    double inner_product(std::vector<double> a, std::vector<double> b);              // inner product of a and b
    double calc_bending_of_two_t(std::vector<double> t1, std::vector<double> t2);    // calculate the bending angle between two tangent

    std::vector<double> Rodrigues_rotation(std::vector<double> v, std::vector<double> k, double theta); // v rotate around k by theta
};
#endif