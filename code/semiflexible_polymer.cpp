#include "semiflexible_polymer.h"
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
// #define PI 3.14159265358979323846

// initialization
semiflexible_polymer::semiflexible_polymer(double L_, Energy_parameter Epar_, double d_theta_, int fixe_bead_0_)
{
    // system related
    beta = 1;
    L = L_;
    fix_bead_0 = fixe_bead_0_;

    // set energy related
    // geometric
    Epar.kappa = Epar_.kappa;
    Epar.f = Epar_.f;
    Epar.g = Epar_.g;
    d_theta = d_theta_;

    // set random number generators
    std::random_device rd;
    std::mt19937 gen_set(rd());
    std::uniform_real_distribution<> rand_uni_set(0, 1);
    std::normal_distribution<> rand_norm_set(0, 1); // normal

    gen = gen_set;
    rand_uni = rand_uni_set;
    rand_norm = rand_norm_set;

    std::cout << "setting system param:" << "L" << L << ", kappa:" << Epar.kappa << ", f:" << Epar.f << ", g:" << Epar.g << "\n";

    // initial configuration
    polymer.resize(L + 1); // L+1 beads in total
    E_sys = 0;
    Tb_sys = 0;

    // initialize the polymer as straight along x axis,
    // such that there is no bending the shear energy
    for (int i = 0; i < L; i++)
    {
        polymer[i].r = {1.0 * i, 0, 0};
        polymer[i].t = {1.0, 0, 0};
        E_sys += -Epar.f * polymer[i].t[0];
    }
    polymer[L].r = {1.0 * L, 0, 0};
    polymer[L].t = {0, 0, 0}; // no tangent for the last bead

    // above initial configuration give 0 bending and shear energy
}

int semiflexible_polymer::update_bead_concerted_rotation(int bead_i)
{
    // bead_i is between 1 and L-1
    bead old_bead;
    bead old_bead_pre; // this is also affected then i is not on the start/end
    double old_E = 0;
    double new_E = 0;
    double old_Tb = 0;
    double new_Tb = 0;
    int satisfy_self_avoiding = 0;

    old_bead = polymer[bead_i];
    old_bead_pre = polymer[bead_i - 1]; // this is also affected then i is not on the start/end
    // bead_i is the middle (i-1)--(i)--(i+1)--
    // what updates:
    // t(i-1), t(i), r(i)
    // Energy change related: since (i-1).t and (i).t changes
    // bending (i-2)--(i-1)-- on (i-1) and  (i)--(i+1)-- on (i+1)
    // streching: no change since (i-1) and (i+1) are fixed
    // shearing (i-1)-- and (i)--

    old_Tb = (bead_i == L - 1) ? 0 : calc_bending_of_two_t(polymer[bead_i].t, polymer[bead_i + 1].t);
    old_Tb += (bead_i == 1) ? 0 : calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i - 2].t);
    // bead_i -2 non existent if bead_i == 1
    // bead_i +1 non existent if bead_i == L - 1

    old_E = 0.5 * Epar.kappa * old_Tb;
    old_E += -Epar.g * polymer[bead_i - 1].r[2] * polymer[bead_i - 1].t[0];
    old_E += -Epar.g * polymer[bead_i].r[2] * polymer[bead_i].t[0];

    // update
    // 1. find the unit vector uni_t from (i-1) to (i+1)
    std::vector<double> uni_t = {polymer[bead_i - 1].t[0] + polymer[bead_i].t[0],
                                 polymer[bead_i - 1].t[1] + polymer[bead_i].t[1],
                                 polymer[bead_i - 1].t[2] + polymer[bead_i].t[2]};

    double uni_t_norm = std::sqrt(uni_t[0] * uni_t[0] + uni_t[1] * uni_t[1] + uni_t[2] * uni_t[2]);
    uni_t[0] /= uni_t_norm;
    uni_t[1] /= uni_t_norm;
    uni_t[2] /= uni_t_norm;

    // return if uni_t is parallel to t(i-1) or t(i)
    if (std::abs(inner_product(uni_t, polymer[bead_i - 1].t) - 1) < 1e-6)
    {
        // std::cout<<"Error: uni_t is parallel to t(i-1)"<<std::endl;
        return 1;
    }

    // 2. rotate both (i-1).t and (i).t about uni_t for random angle theta in [-d_theta,d_theta]
    double theta = d_theta * (2 * rand_uni(gen) - 1);
    polymer[bead_i - 1].t = Rodrigues_rotation(polymer[bead_i - 1].t, uni_t, theta);
    polymer[bead_i].t = Rodrigues_rotation(polymer[bead_i].t, uni_t, theta);
    polymer[bead_i].r = {polymer[bead_i - 1].r[0] + polymer[bead_i - 1].t[0],
                         polymer[bead_i - 1].r[1] + polymer[bead_i - 1].t[1],
                         polymer[bead_i - 1].r[2] + polymer[bead_i - 1].t[2]};
    // check for r(i+1)

    /*
    if (bead_i != L-1 && ([bead_i].r[0] + polymer[bead_i].t[0] != polymer[bead_i + 1].r[0] || polymer[bead_i].r[1] + polymer[bead_i].t[1] != polymer[bead_i + 1].r[1] || polymer[bead_i].r[2] + polymer[bead_i].t[2] != polymer[bead_i + 1].r[2]))
    {
        std::cout << "Error: r(i+1) not match after update bead_i" << std::endl;
        return 0;
    }
    */

    // new energy
    new_Tb = (bead_i == L - 1) ? 0 : calc_bending_of_two_t(polymer[bead_i].t, polymer[bead_i + 1].t);
    new_Tb += (bead_i == 1) ? 0 : calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i - 2].t);
    new_E = 0.5 * Epar.kappa * new_Tb;
    new_E += -Epar.g * polymer[bead_i - 1].r[2] * polymer[bead_i - 1].t[0];
    new_E += -Epar.g * polymer[bead_i].r[2] * polymer[bead_i].t[0];

    satisfy_self_avoiding = satisfy_self_avoiding_condition(bead_i);

    // Metropolis
    if (satisfy_self_avoiding && rand_uni(gen) < std::exp(-beta * (new_E - old_E)))
    {
        E_sys += new_E - old_E;
        Tb_sys += new_Tb - old_Tb;
        return 1;
    }
    else
    {
        polymer[bead_i] = old_bead;
        polymer[bead_i - 1] = old_bead_pre;
        return 0;
    }
}

int semiflexible_polymer::update_bead_tangent_rotation(int bead_i)
{
    std::vector<bead> old_beads{}; // starting with bead_i
    double old_E = 0;
    double new_E = 0;
    double old_Tb = 0;
    double new_Tb = 0;
    for (int j = bead_i; j < polymer.size(); j++)
    {
        old_beads.push_back(polymer[j]); // (i),(i+1)... up to (L-1),(L)
    }

    // for any bead_i with a t, which is [0,L-1]
    // everything j > bead_i us updated, since t_i rotate
    // what updates:
    // t(i), r(i+1), t(i+1), ..... r(L-1), t(L-1), r(L)
    // Energ change related
    // bending (i-1)--(i)--(i+1)-- on (i) if i>1 else, no bending changed
    // streching r(L) - r(i)
    // shearing all j>bead_i

    old_Tb = (bead_i == 0) ? 0 : calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i].t);
    old_E = 0.5 * Epar.kappa * old_Tb;
    old_E += -Epar.f * (polymer[L].r[0] - polymer[bead_i].r[0]);
    for (int j = bead_i; j < polymer.size(); j++)
    {
        old_E += -Epar.g * polymer[j].r[2] * polymer[j].t[0];
    }

    // update
    // 1. find random unit uni_v vector not parallel to t
    std::vector<double> uni_v = rand_uni_vec();
    while (std::abs(inner_product(uni_v, polymer[bead_i].t)) < 1e-6)
    { // uni_v must be not parallel to t
        uni_v = rand_uni_vec();
    }
    // 2. transform uni_v to vector uni_w perpendicular to t
    std::vector<double> uni_w = cross_product(polymer[bead_i].t, uni_v);
    double uni_w_norm = std::sqrt(uni_w[0] * uni_w[0] + uni_w[1] * uni_w[1] + uni_w[2] * uni_w[2]);
    uni_w[0] /= uni_w_norm;
    uni_w[1] /= uni_w_norm;
    uni_w[2] /= uni_w_norm;

    // also need uni_v perpendicular to both w and t
    uni_v = cross_product(polymer[bead_i].t, uni_w);
    double uni_v_norm = std::sqrt(uni_v[0] * uni_v[0] + uni_v[1] * uni_v[1] + uni_v[2] * uni_v[2]);
    uni_v[0] /= uni_v_norm;
    uni_v[1] /= uni_v_norm;
    uni_v[2] /= uni_v_norm;

    // 3. rotate every r(ij) towards uni_w for random angle theta such that cos(theta) uniform in [cos d_theta,1]
    // due to the length of rotated polymer, here use d_theta/(L-bead_i), to adjust acceptance rate
    double cos_theta = 1 - rand_uni(gen) * (1 - std::cos(1.0 * d_theta / (1 + (Epar.f + Epar.g) * (L - bead_i))));
    double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
    std::vector<double> r_ij(3, 0);
    double r_ij_norm = 0;
    double r_ij_norm_new = 0;
    double r_ij_v = 0;
    double r_ij_w = 0;
    double r_ij_t = 0;
    // note: txw = v
    // check othorgonal of v,w,t
    if (std::abs(inner_product(uni_v, uni_w)) > 1e-6 || std::abs(inner_product(uni_v, polymer[bead_i].t)) > 1e-6 || std::abs(inner_product(uni_w, polymer[bead_i].t)) > 1e-6)
    {
        std::cout << "Error: uni_v, uni_w, t not orthogonal" << std::endl;
        return 0;
    }

    for (int j = bead_i + 1; j < polymer.size(); j++)
    {
        // rorate r_ij towards uni_w by theta
        r_ij = {polymer[j].r[0] - polymer[bead_i].r[0], polymer[j].r[1] - polymer[bead_i].r[1], polymer[j].r[2] - polymer[bead_i].r[2]};
        r_ij_norm = std::sqrt(r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2]);
        r_ij_v = inner_product(r_ij, uni_v);
        r_ij_w = inner_product(r_ij, uni_w);
        r_ij_t = inner_product(r_ij, polymer[bead_i].t);

        r_ij = {(cos_theta * r_ij_t - sin_theta * r_ij_w) * polymer[bead_i].t[0] + (sin_theta * r_ij_t + cos_theta * r_ij_w) * uni_w[0] + r_ij_v * uni_v[0],
                (cos_theta * r_ij_t - sin_theta * r_ij_w) * polymer[bead_i].t[1] + (sin_theta * r_ij_t + cos_theta * r_ij_w) * uni_w[1] + r_ij_v * uni_v[1],
                (cos_theta * r_ij_t - sin_theta * r_ij_w) * polymer[bead_i].t[2] + (sin_theta * r_ij_t + cos_theta * r_ij_w) * uni_w[2] + r_ij_v * uni_v[2]};

        // renormalize r_ij
        r_ij_norm_new = std::sqrt(r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2]);
        r_ij[0] *= r_ij_norm / r_ij_norm_new;
        r_ij[1] *= r_ij_norm / r_ij_norm_new;
        r_ij[2] *= r_ij_norm / r_ij_norm_new;

        polymer[j].r = {polymer[bead_i].r[0] + r_ij[0], polymer[bead_i].r[1] + r_ij[1], polymer[bead_i].r[2] + r_ij[2]};

        // check is r_ij_norm changed
        if (std::abs(r_ij_norm - std::sqrt(r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2])) > 1e-6)
        {
            std::cout << "Error: r_ij_norm changed" << std::endl;
            // return 0;
        }
    }

    // update all t accodingly
    for (int j = bead_i; j < polymer.size() - 1; j++)
    {
        polymer[j].t = {polymer[j + 1].r[0] - polymer[j].r[0], polymer[j + 1].r[1] - polymer[j].r[1], polymer[j + 1].r[2] - polymer[j].r[2]};
        // check if t is normalized
        if (std::abs(polymer[j].t[0] * polymer[j].t[0] + polymer[j].t[1] * polymer[j].t[1] + polymer[j].t[2] * polymer[j].t[2] - 1) > 1e-6)
        {
            std::cout << "Error: t not normalized" << std::endl;
        }
    }

    // calculate updated energy
    new_Tb = (bead_i == 0) ? 0 : calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i].t);
    new_E = 0.5 * Epar.kappa * new_Tb;
    new_E += -Epar.f * (polymer[L].r[0] - polymer[bead_i].r[0]);
    for (int j = bead_i; j < polymer.size(); j++)
    {
        new_E += -Epar.g * polymer[j].r[2] * polymer[j].t[0];
    }

    if (satisfy_self_avoiding_condition_by_group(bead_i) && rand_uni(gen) < std::exp(-beta * (new_E - old_E)))
    {
        E_sys += new_E - old_E;
        Tb_sys += new_Tb - old_Tb;
        return 1;
    }
    else
    {
        for (int j = bead_i; j < polymer.size(); j++)
        {
            polymer[j] = old_beads[j - bead_i];
        }
        return 0;
    }
}

int semiflexible_polymer::satisfy_self_avoiding_condition(int bead_i)
{
    double rij_norm2 = 0;
    std::vector<double> bead_r(3, 0);
    if (bead_i == L)
    {
        bead_r[0] = polymer[bead_i - 1].r[0] + polymer[bead_i - 1].t[0];
        bead_r[1] = polymer[bead_i - 1].r[1] + polymer[bead_i - 1].t[1];
        bead_r[2] = polymer[bead_i - 1].r[2] + polymer[bead_i - 1].t[2];
    }
    else
    {
        bead_r = polymer[bead_i].r;
    }
    // check distance between bead_i and all other beads
    for (int j = 0; j < L; j++)
    {
        if (j == bead_i - 1 || j == bead_i || j == bead_i + 1)
        {
            // neighbouring beads and self not need to check
            continue;
        }
        rij_norm2 = 0;
        for (int k = 0; k < 3; k++)
        {
            rij_norm2 += std::pow(bead_r[k] - polymer[j].r[k], 2);
        }
        // comparing the distance between i's next segment's center and j
        if (rij_norm2 < 1)
        {
            return 0;
        }
    }
    return 1;
}

int semiflexible_polymer::satisfy_self_avoiding_condition_by_group(int bead_i)
{
    // check overlapping between j<bead_i and j>bead_i
    if (bead_i == 0 || bead_i == L)
    {
        return 1;
    }
    double rjk_norm2 = 0;
    for (int j = 0; j < bead_i; j++)
    {
        for (int k = bead_i + 1; k < polymer.size(); k++)
        {
            rjk_norm2 = 0;
            for (int l = 0; l < 3; l++)
            {
                rjk_norm2 += std::pow(polymer[j].r[l] - polymer[k].r[l], 2);
            }
            if (rjk_norm2 < 1)
            {
                return 0;
            }
        }
    }
    return 1;
}

void semiflexible_polymer::save_polymer_to_file(std::string filename)
{
    std::ofstream f(filename);
    if (f.is_open())
    {
        f << "r[0],r[1],r[2],t[0],t[1],t[2]";

        for (int i = 0; i < size(polymer); i++)
        {
            f << "\n"
              << polymer[i].r[0] << "," << polymer[i].r[1] << "," << polymer[i].r[2]
              << "," << polymer[i].t[0] << "," << polymer[i].t[1] << "," << polymer[i].t[2];
        }
    }
    f.close();

} // end save_polymer_to_file

void semiflexible_polymer::save_observable_to_file(std::string filename, std::vector<observable> obs_ensemble, bool save_detail)
{
    std::ofstream f(filename);
    if (f.is_open())
    {
        int number_of_polymer = obs_ensemble.size();
        if (save_detail)
        {
            f << "beta=" << beta << "\n";
            f << "f=" << Epar.f << "\n";
            f << "g=" << Epar.g << "\n";
            f << "E(energy),Tb(total bending), X(end to end X distance),  Y(end to end Y distance), Z(end to end Z distance), R(end to end distance)";
            for (int i = 0; i < number_of_polymer; i++)
            {
                f << "\n"
                  << obs_ensemble[i].E << "," << obs_ensemble[i].Tb
                  << "," << obs_ensemble[i].X << "," << obs_ensemble[i].Y
                  << "," << obs_ensemble[i].Z << "," << obs_ensemble[i].R;
            }
        }
        else
        {
            // only save stats
            // find stats
            // calculate average and standard deviation of E
            double avg_E = 0.0;
            double std_E = 0.0;
            double avg_Tb = 0.0;
            double std_Tb = 0.0;
            double avg_X = 0.0;
            double std_X = 0.0;
            double avg_Y = 0.0;
            double std_Y = 0.0;
            double avg_Z = 0.0;
            double std_Z = 0.0;
            double avg_R = 0.0;
            double std_R = 0.0;
            std::vector<double> avg_Sq(obs_ensemble[0].Sq.size(), 0.0);
            std::vector<double> std_Sq(obs_ensemble[0].Sq.size(), 0.0);
            std::vector<double> avg_tts(obs_ensemble[0].tts.size(), 0.0);
            std::vector<double> std_tts(obs_ensemble[0].tts.size(), 0.0);

            // get the average
            for (int i = 0; i < obs_ensemble.size(); i++)
            {
                avg_E += obs_ensemble[i].E;
                avg_Tb += obs_ensemble[i].Tb;
                avg_X += obs_ensemble[i].X;
                avg_Y += obs_ensemble[i].Y;
                avg_Z += obs_ensemble[i].Z;
                avg_R += obs_ensemble[i].R;
                for (int j = 0; j < obs_ensemble[i].Sq.size(); j++)
                {
                    avg_Sq[j] += obs_ensemble[i].Sq[j];
                    avg_tts[j] += obs_ensemble[i].tts[j];
                }
            }
            avg_E /= obs_ensemble.size();
            avg_Tb /= obs_ensemble.size();
            avg_X /= obs_ensemble.size();
            avg_Y /= obs_ensemble.size();
            avg_Z /= obs_ensemble.size();
            avg_R /= obs_ensemble.size();

            for (int j = 0; j < obs_ensemble[0].Sq.size(); j++)
            {
                avg_Sq[j] /= obs_ensemble.size();
                avg_tts[j] /= obs_ensemble.size();
            }

            // get std

            for (int i = 0; i < obs_ensemble.size(); i++)
            {
                std_E += (obs_ensemble[i].E - avg_E) * (obs_ensemble[i].E - avg_E);
                std_Tb += (obs_ensemble[i].Tb - avg_Tb) * (obs_ensemble[i].Tb - avg_Tb);
                std_X += (obs_ensemble[i].X - avg_X) * (obs_ensemble[i].X - avg_X);
                std_Y += (obs_ensemble[i].Y - avg_Y) * (obs_ensemble[i].Y - avg_Y);
                std_Z += (obs_ensemble[i].Z - avg_Z) * (obs_ensemble[i].Z - avg_Z);
                std_R += (obs_ensemble[i].R - avg_R) * (obs_ensemble[i].R - avg_R);

                for (int j = 0; j < obs_ensemble[i].Sq.size(); j++)
                {
                    std_Sq[j] += (obs_ensemble[i].Sq[j] - avg_Sq[j]) * (obs_ensemble[i].Sq[j] - avg_Sq[j]);
                    std_tts[j] += (obs_ensemble[i].tts[j] - avg_tts[j]) * (obs_ensemble[i].tts[j] - avg_tts[j]);
                }
            }
            std_E = std::sqrt(std_E / obs_ensemble.size());
            std_Tb = std::sqrt(std_Tb / obs_ensemble.size());
            std_X = std::sqrt(std_X / obs_ensemble.size());
            std_Y = std::sqrt(std_Y / obs_ensemble.size());
            std_Z = std::sqrt(std_Z / obs_ensemble.size());
            std_R = std::sqrt(std_R / obs_ensemble.size());
            for (int j = 0; j < obs_ensemble[0].Sq.size(); j++)
            {
                std_Sq[j] = std::sqrt(std_Sq[j] / obs_ensemble.size());
                std_tts[j] = std::sqrt(std_tts[j] / obs_ensemble.size());
            }

            // write parameters and stats to the file
            f << "stats,L,kappa,f,g,E,Tb,X,Y,Z,R,Sq";
            f << "\nmean," << L << "," << Epar.kappa << "," << Epar.f << "," << Epar.g << "," << avg_E << "," << avg_Tb
              << "," << avg_X << "," << avg_Y << "," << avg_Z << "," << avg_R;
            for (int j = 0; j < avg_Sq.size(); j++)
            {
                f << "," << avg_Sq[j];
            }
            f << "\nstd/sqrt(number of polymer),NA,NA,NA,NA," << std_E / std::sqrt(obs_ensemble.size()) << "," << std_Tb / std::sqrt(obs_ensemble.size())
              << "," << std_X / std::sqrt(obs_ensemble.size()) << "," << std_Y / std::sqrt(obs_ensemble.size())
              << "," << std_Z / std::sqrt(obs_ensemble.size()) << "," << std_R / std::sqrt(obs_ensemble.size());
            for (int j = 0; j < std_Sq.size(); j++)
            {
                f << "," << std_Sq[j] / std::sqrt(obs_ensemble.size());
            }

            f << "\n qB,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < obs_ensemble[0].qB.size(); j++)
            {
                f << "," << obs_ensemble[0].qB[j];
            }

            f << "\n tts,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < avg_tts.size(); j++)
            {
                f << "," << avg_tts[j];
            }
            f << "\nstd_tts/sqrt(number of polymer),NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < std_tts.size(); j++)
            {
                f << "," << std_tts[j] / std::sqrt(obs_ensemble.size());
            }
            f << "\n r(for tts),NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < obs_ensemble[0].spB.size(); j++)
            {
                f << "," << obs_ensemble[0].spB[j];
            }
        }
    }
    f.close();
}

observable semiflexible_polymer::measure_observable(int bin_num)
{
    // measure observable
    observable obs;
    obs.E = E_sys;
    obs.Tb = Tb_sys;

    std::vector<double> R = {0, 0, 0};
    R[0] = polymer[L - 1].r[0] + polymer[L - 1].t[0] - polymer[0].r[0];
    R[1] = polymer[L - 1].r[1] + polymer[L - 1].t[1] - polymer[0].r[1];
    R[2] = polymer[L - 1].r[2] + polymer[L - 1].t[2] - polymer[0].r[2];

    obs.R = std::sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
    obs.X = R[0];
    obs.Y = R[1];
    obs.Z = R[2];

    double qB_i = 1e-4; // 0.2*M_PI/L; //0.1/L; ;
    double qB_f = 1e0;  // M_PI;//100.0/L; //M_PI;
    obs.qB = std::vector<double>(bin_num, 0);
    for (int k = 0; k < bin_num; k++)
    {
        obs.qB[k] = qB_i * std::pow(qB_f / qB_i, 1.0 * k / (bin_num - 1)); // uniform in log scale
    }
    obs.Sq = calc_structure_factor(obs.qB);


    obs.spB = std::vector<double>(int(L/4), 0); // measure upto L/5 length
    for (int k = 0; k < int(L/4); k++)
    {
        obs.spB[k] = k; // 0, 1, 2, ... bin_num-1
    }
    obs.tts = calc_tangent_pair_correlation(obs.spB);

    return obs;
}

std::vector<double> semiflexible_polymer::calc_structure_factor(std::vector<double> qB)
{
    // measure structure factor
    int bin;
    double q;
    int bin_num = qB.size();

    std::vector<std::vector<double>> R_all{}; // all scattering point's R[axis] value, including all scattering points in each segment

    for (int i = 0; i < size(polymer); i++)
    {
        R_all.push_back(polymer[i].r);
    }
    int N = R_all.size();                      // total number of scattering points
    double r;                                  // for storating distance between two scattering points
    std::vector<double> SqB(bin_num, 1.0 / N); // initialize with 1 due to self overlaping term (see S.H. Chen 1986 eq 18)
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            r = std::sqrt(std::pow(R_all[i][0] - R_all[j][0], 2) + std::pow(R_all[i][1] - R_all[j][1], 2) + std::pow(R_all[i][2] - R_all[j][2], 2));
            // calculate S_q
            for (int k = 0; k < bin_num; k++)
            {
                q = qB[k];                                         // B = 1
                SqB[k] += 2.0 / N / N * std::sin(q * r) / (q * r); // normalize to 1 at q=0
            }
        }
    }
    return SqB;
}

std::vector<double> semiflexible_polymer::calc_tangent_pair_correlation(std::vector<double> spB)
{
    int bin_num = spB.size();
    std::vector<double> tts(bin_num, 0);
    tts[0] = 1.0;                            // self correlation
    std::vector<int> pair_count(bin_num, 0); // count number of t-t pair added to a s
    for (int i = 0; i < L - 1; i++)
    {
        for (int j = i + 1; j < L; j++)
        {
            tts[j - i] += inner_product(polymer[i].t, polymer[j].t);
            pair_count[j - i]++;
        }
    }
    for (int i = 1; i < bin_num; i++)
    {
        tts[i] /= pair_count[i];
    }
    return tts;
}

double semiflexible_polymer::calc_radius_of_gyration()
{
    // calculate radius of gyration
    double Rg = 0;
    // find center of mass
    double x_c = 0, y_c = 0, z_c = 0;
    for (int i = 0; i < L; i++)
    {
        x_c = x_c + polymer[i].r[0];
        y_c = y_c + polymer[i].r[1];
        z_c = z_c + polymer[i].r[2];
    }
    x_c = x_c / L;
    y_c = y_c / L;
    z_c = z_c / L;

    // calculate Rg
    for (int i = 0; i < L; i++)
    {
        Rg += std::pow(polymer[i].r[0] - x_c, 2) + std::pow(polymer[i].r[1] - y_c, 2) + std::pow(polymer[i].r[2] - z_c, 2);
    }
    Rg = std::sqrt(Rg / L);
    return Rg;
}

#pragma region : useful tools
std::vector<double> semiflexible_polymer::calc_rod_structure_factor(std::vector<double> qB)
{
    int L = size(polymer);
    // measure structure factor
    std::vector<double> SqB(qB.size(), 1.0 / L); // initialize with 1 due to self overlaping term (see S.H. Chen 1986 eq 18)
    for (int i = 0; i < L - 1; i++)
    {
        for (int j = i + 1; j < L; j++)
        {
            for (int k = 0; k < qB.size(); k++)
            {
                SqB[k] += 2.0 / L / L * std::sin(qB[k] * (j - i)) / (qB[k] * (j - i)); // normalize to 1 at q=0
            }
        }
    }
    return SqB;
}

std::vector<double> semiflexible_polymer::cross_product(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> c(3, 0);
    for (int j = 0; j < 3; j++)
    {
        c[j] = a[(j + 1) % 3] * b[(j + 2) % 3] - a[(j + 2) % 3] * b[(j + 1) % 3];
    }
    return c;
}

double semiflexible_polymer::inner_product(std::vector<double> a, std::vector<double> b)
{
    double c = 0;
    for (int j = 0; j < 3; j++)
    {
        c += a[j] * b[j];
    }
    return c;
}

std::vector<double> semiflexible_polymer::Rodrigues_rotation(std::vector<double> v, std::vector<double> k, double theta)
{
    // v rotate around k by theta
    std::vector<double> v_rot(3, 0);
    // 1. find (k x v), k cross v and kv, k.v
    std::vector<double> kxv = cross_product(k, v);
    double kv = inner_product(k, v);
    // 2. do the calculation
    for (int j = 0; j < 3; j++)
    {
        v_rot[j] = v[j] * std::cos(theta) + kxv[j] * std::sin(theta) + k[j] * kv * (1 - std::cos(theta));
    }
    // 3. renormalize v_rot
    double v_rot_norm = std::sqrt(v_rot[0] * v_rot[0] + v_rot[1] * v_rot[1] + v_rot[2] * v_rot[2]);
    v_rot[0] /= v_rot_norm;
    v_rot[1] /= v_rot_norm;
    v_rot[2] /= v_rot_norm;

    if (std::abs(v_rot[0] * v_rot[0] + v_rot[1] * v_rot[1] + v_rot[2] * v_rot[2] - 1) > 1e-6)
    {
        printf("v_rot^2=%f\n", v_rot[0] * v_rot[0] + v_rot[1] * v_rot[1] + v_rot[2] * v_rot[2]);
        throw std::runtime_error("Error: v_rot^2 is not 1.");
    }

    return v_rot;
}

std::vector<double> semiflexible_polymer::rand_uni_vec()
{
    // 1. generat random point in unit ball
    double x, y, z;
    do
    {
        x = 2 * rand_uni(gen) - 1;
        y = 2 * rand_uni(gen) - 1;
        z = 2 * rand_uni(gen) - 1;
    } while (x * x + y * y + z * z > 1);
    // 2. normalize
    double norm = std::sqrt(x * x + y * y + z * z);
    x /= norm;
    y /= norm;
    z /= norm;
    return {x, y, z};
}

double semiflexible_polymer::calc_bending_of_two_t(std::vector<double> t1, std::vector<double> t2)
{
    double ans = 0;
    ans += (t1[0] - t2[0]) * (t1[0] - t2[0]);
    ans += (t1[1] - t2[1]) * (t1[1] - t2[1]);
    ans += (t1[2] - t2[2]) * (t1[2] - t2[2]);
    return ans;
}

void semiflexible_polymer::run_simultion(int therm_sweep, int MC_sweeps, int step_per_sweep, std::string folder, std::string finfo, int bin_num)
{
    std::vector<observable> obs_ensemble;

    double conrot_acceptance_rate = 0;
    double tanrot_acceptance_rate = 0;
    int bead_i;

    beta = 0; // randomization
    for (int i = 0; i < therm_sweep; i++)
    {
        std::cout << "randomization(beta=0) sweep " << i << " out of " << therm_sweep << " (" << (i * 100) / therm_sweep << "%)\r";
        for (int j = 0; j < step_per_sweep; j++)
        {
            bead_i = 1 + int(rand_uni(gen) * (L - 2));
            update_bead_concerted_rotation(bead_i);
        }
        bead_i = int(rand_uni(gen) * L);
        update_bead_tangent_rotation(bead_i);
    }

    std::cout << "\n";
    // thermalization using tempering
    for (int i = 0; i < therm_sweep; i++)
    {
        beta = (i + 1) * 1.0 / therm_sweep; // tempering
        std::cout << "thermalization/tempering (beta -> 1) sweep " << i << " out of " << therm_sweep << " (" << (i * 100) / therm_sweep << "%)\r";
        for (int j = 0; j < step_per_sweep; j++)
        {
            bead_i = 1 + int(rand_uni(gen) * (L - 2));
            update_bead_concerted_rotation(bead_i);
        }
        bead_i = int(rand_uni(gen) * L);
        update_bead_tangent_rotation(bead_i);
    }

    // data collection run
    // TODO: consider implementing irreversible algorithem like event-chain here if too slow
    std::cout << "\n";
    if (beta != 1)
    {
        std::cout << "Error: beta is not 1 at the end of thermalization\n";
    }
    for (int i = 0; i < MC_sweeps; i++)
    {
        std::cout << "MC sweep " << i << " out of " << MC_sweeps << " (" << (i * 100) / MC_sweeps << "%)\r";
        for (int j = 0; j < step_per_sweep; j++)
        {
            bead_i = 1 + int(rand_uni(gen) * (L - 2));
            if (fix_bead_0 && bead_i == 0)
            {
                std::cout << "Error: bead_i is 0 when fixed_bead_0" << std::endl;
            }
            if (bead_i == L + 1)
            {
                std::cout << "bead_i == L + 1????" << std::endl;
            }
            // std::cout<<"bead_i:"<<bead_i<<std::endl;
            // if(bead_i == L){
            // std::cout<<"updating the end"<<std::endl;
            //}
            conrot_acceptance_rate += update_bead_concerted_rotation(bead_i);
        }
        bead_i = int(rand_uni(gen) * L);
        tanrot_acceptance_rate += update_bead_tangent_rotation(bead_i);
        obs_ensemble.push_back(measure_observable(bin_num));
    }
    std::cout << "\n";
    std::cout << "conrot_acceptance rate:" << conrot_acceptance_rate / MC_sweeps / step_per_sweep << std::endl;
    std::cout << "tanrot_acceptance rate:" << tanrot_acceptance_rate / MC_sweeps << std::endl;

    save_polymer_to_file(folder + "/config_" + finfo + ".csv");
    save_observable_to_file(folder + "/obs_MC_" + finfo + ".csv", obs_ensemble, true);
    save_observable_to_file(folder + "/obs_" + finfo + ".csv", obs_ensemble, false);
}