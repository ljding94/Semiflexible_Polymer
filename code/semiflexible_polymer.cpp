#include "semiflexible_polymer.h"
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <filesystem>
// #define PI 3.14159265358979323846

// initialization
semiflexible_polymer::semiflexible_polymer(double L_, Energy_parameter Epar_, bool random_Epar)
{
    // system related
    beta = 1;
    L = L_;

    // set energy related
    // geometric

    // set random number generators
    std::random_device rd;
    std::mt19937 gen_set(rd());
    std::uniform_real_distribution<> rand_uni_set(0, 1);
    std::normal_distribution<> rand_norm_set(0, 1); // normal

    gen = gen_set;
    rand_uni = rand_uni_set;
    rand_norm = rand_norm_set;

    if (random_Epar)
    {
        Epar.kappa = 2 + 18 * rand_uni(gen);
        Epar.f = 0.5 * rand_uni(gen);
        Epar.g = 2.0 * rand_uni(gen)/L;
    }
    else
    {
        Epar.kappa = Epar_.kappa;
        Epar.f = Epar_.f;
        Epar.g = Epar_.g;
    }

    d_theta = M_PI * 2.0 / 3 / (1.0 + std::sqrt(Epar.kappa));
    std::cout << "setting system param:" << "L" << L << ", kappa:" << Epar.kappa << ", f:" << Epar.f << ", g:" << Epar.g << "\n";

    // initial configuration
    polymer.resize(L + 1); // L+1 beads in total
    E_sys = 0;
    Tb_sys = 0;

    // initialize the polymer as straight along x axis,
    // such that there is no bending the shear energy
    // int mid_point = int(L / 2);
    for (int i = 0; i < L; i++)
    {
        polymer[i].r = {0, 0, 1.0 * i};
        polymer[i].t = {0, 0, 1.0};
        E_sys += -Epar.f * polymer[i].t[0];
    }
    polymer[L].r = {0, 0, 1.0 * L};
    polymer[L].t = {0, 0, 0}; // no tangent for the last bead

    // above initial configuration give 0 bending and shear energy
}

int semiflexible_polymer::update_bead_crankshaft(int bead_i, int bead_j)
{
    if (bead_i > bead_j)
    {
        // swap, make sure i is before j
        int bead = bead_i;
        bead_i = bead_j;
        bead_j = bead;
    }
    bead_j = std::min(bead_j, L);
    if (bead_j <= bead_i + 2)
    {
        return 0;
    }
    std::vector<bead> old_beads{}; // [bead_i, bead_j]
    double old_E = 0;
    double new_E = 0;
    double old_Tb = 0;
    double new_Tb = 0;

    for (int j = bead_i; j <= bead_j; j++)
    {
        old_beads.push_back(polymer[j]); // (i),(i+1)... up to --(j),
    }
    // what updates:
    // t(i), r(i+1), t(i+1),..... r(j-1), t(j-1)
    // Energy change related:
    // bending (i-1)--(i)-- on (i) and  (j-1)--(j)-- on (i)
    // streching: no change since (0) and (L) are fixed
    // shearing (i) up to (j-1)

    old_Tb = (bead_i == 0) ? 0 : calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i].t);
    old_Tb += (bead_j == L) ? 0 : calc_bending_of_two_t(polymer[bead_j - 1].t, polymer[bead_j].t);
    // bead_L has no t

    old_E = 0.5 * Epar.kappa * old_Tb;
    for (int j = bead_i; j < bead_j; j++)
    {
        old_E += -Epar.g * polymer[j].r[2] * polymer[j].t[0];
    }

    // update
    // 1. find the unit vector uni_t from (i) to (j)
    std::vector<double> uni_t = {polymer[bead_j].r[0] - polymer[bead_i].r[0],
                                 polymer[bead_j].r[1] - polymer[bead_i].r[1],
                                 polymer[bead_j].r[2] - polymer[bead_i].r[2]};

    double uni_t_norm = std::sqrt(uni_t[0] * uni_t[0] + uni_t[1] * uni_t[1] + uni_t[2] * uni_t[2]);
    uni_t[0] /= uni_t_norm;
    uni_t[1] /= uni_t_norm;
    uni_t[2] /= uni_t_norm;

    // 2. rotate everything between (i) and (j) about uni_t for random angle theta in renormalized [-d_theta,d_theta]
    double theta = d_theta / (1 + (Epar.f + Epar.g) * (bead_j - bead_i)) * (2 * rand_uni(gen) - 1);
    // TODO: consider using simpler algorithm than Rodrigues_rotation

    std::vector<double> rik{0, 0, 0};
    double rik_norm;
    for (int k = bead_i; k < bead_j; k++)
    {
        // rik = {polymer[k].r[0] - polymer[bead_i].r[0], polymer[k].r[1] - polymer[bead_i].r[1], polymer[k].r[2] - polymer[bead_i].r[2]};
        // rik_norm = std::sqrt(rik[0] * rik[0] + rik[1] * rik[1] + rik[2] * rik[2]);
        polymer[k].t = Rodrigues_rotation(polymer[k].t, uni_t, theta); // this returned renomalized result
    }
    // reconstruct r
    for (int k = bead_i + 1; k < bead_j; k++)
    {
        polymer[k].r = {polymer[k - 1].r[0] + polymer[k - 1].t[0],
                        polymer[k - 1].r[1] + polymer[k - 1].t[1],
                        polymer[k - 1].r[2] + polymer[k - 1].t[2]};
    }
    if (std::abs(polymer[bead_j - 1].r[0] + polymer[bead_j - 1].t[0] - polymer[bead_j].r[0]) > 1e-6 ||
        std::abs(polymer[bead_j - 1].r[1] + polymer[bead_j - 1].t[1] - polymer[bead_j].r[1]) > 1e-6 ||
        std::abs(polymer[bead_j - 1].r[2] + polymer[bead_j - 1].t[2] - polymer[bead_j].r[2]) > 1e-6)
    {
        std::cout << "Error: r(j) not match after update bead_i" << std::endl;
        std::cout << "0:" << polymer[bead_j - 1].r[0] + polymer[bead_j - 1].t[0] << " " << polymer[bead_j].r[0] << std::endl;
        std::cout << "1:" << polymer[bead_j - 1].r[1] + polymer[bead_j - 1].t[1] << " " << polymer[bead_j].r[1] << std::endl;
        std::cout << "2:" << polymer[bead_j - 1].r[2] + polymer[bead_j - 1].t[2] << " " << polymer[bead_j].r[2] << std::endl;
        return 0;
    }
    // new energy
    new_Tb = (bead_i == 0) ? 0 : calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i].t);
    new_Tb += (bead_j == L) ? 0 : calc_bending_of_two_t(polymer[bead_j - 1].t, polymer[bead_j].t);
    // bead_L has no t
    new_E = 0.5 * Epar.kappa * new_Tb;
    for (int j = bead_i; j < bead_j; j++)
    {
        new_E += -Epar.g * polymer[j].r[2] * polymer[j].t[0];
    }

    // Metropolis
    if (satisfy_self_avoiding_condition_by_group(bead_i, bead_j) && rand_uni(gen) < std::exp(-beta * (new_E - old_E)))
    {
        E_sys += new_E - old_E;
        Tb_sys += new_Tb - old_Tb;
        return 1;
    }
    else
    {
        for (int j = bead_i; j <= bead_j; j++)
        {
            polymer[j] = old_beads[j - bead_i];
        }

        return 0;
    }
}

int semiflexible_polymer::update_bead_pivot_right(int bead_i)
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
    double d_theta_in_use = (bead_i == 0 && Epar.f == 0 && Epar.g == 0) ? d_theta : d_theta / (1 + (Epar.f + Epar.g) * (L - bead_i));
    double cos_theta = 1 - rand_uni(gen) * (1 - std::cos(1.0 * d_theta_in_use));
    double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
    std::vector<double> r_ij(3, 0);
    double r_ij_norm = 0;
    double r_ij_norm_new = 0;
    double r_ij_v = 0;
    double r_ij_w = 0;
    double r_ij_t = 0;

    double tj_v = 0;
    double tj_w = 0;
    double tj_t = 0;
    double tj_norm = 0;
    // note: txw = v
    // check othorgonal of v,w,t
    if (std::abs(inner_product(uni_v, uni_w)) > 1e-6 || std::abs(inner_product(uni_v, polymer[bead_i].t)) > 1e-6 || std::abs(inner_product(uni_w, polymer[bead_i].t)) > 1e-6)
    {
        std::cout << "Error: uni_v, uni_w, t not orthogonal" << std::endl;
        std::cout << inner_product(uni_v, uni_w) << ":" << inner_product(uni_v, polymer[bead_i].t) << ":" << inner_product(uni_w, polymer[bead_i].t) << std::endl;
        return 0;
    }
    std::vector<double> ti_old = polymer[bead_i].t;
    for (int j = bead_i; j < polymer.size() - 1; j++) // polymer[L] has no t
    {
        // rotate tj towards uni_w by theta
        tj_v = inner_product(polymer[j].t, uni_v);
        tj_w = inner_product(polymer[j].t, uni_w);
        tj_t = inner_product(polymer[j].t, ti_old);
        polymer[j].t = {(cos_theta * tj_t - sin_theta * tj_w) * ti_old[0] + (sin_theta * tj_t + cos_theta * tj_w) * uni_w[0] + tj_v * uni_v[0],
                        (cos_theta * tj_t - sin_theta * tj_w) * ti_old[1] + (sin_theta * tj_t + cos_theta * tj_w) * uni_w[1] + tj_v * uni_v[1],
                        (cos_theta * tj_t - sin_theta * tj_w) * ti_old[2] + (sin_theta * tj_t + cos_theta * tj_w) * uni_w[2] + tj_v * uni_v[2]};
        // renormalize t
        tj_norm = std::sqrt(polymer[j].t[0] * polymer[j].t[0] + polymer[j].t[1] * polymer[j].t[1] + polymer[j].t[2] * polymer[j].t[2]);
        polymer[j].t[0] /= tj_norm;
        polymer[j].t[1] /= tj_norm;
        polymer[j].t[2] /= tj_norm;
    }

    // reconstruct r
    for (int k = bead_i + 1; k < polymer.size(); k++)
    {
        polymer[k].r = {polymer[k - 1].r[0] + polymer[k - 1].t[0],
                        polymer[k - 1].r[1] + polymer[k - 1].t[1],
                        polymer[k - 1].r[2] + polymer[k - 1].t[2]};
    }

    // calculate updated energy
    new_Tb = (bead_i == 0) ? 0 : calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i].t);
    new_E = 0.5 * Epar.kappa * new_Tb;
    new_E += -Epar.f * (polymer[L].r[0] - polymer[bead_i].r[0]);
    for (int j = bead_i; j < polymer.size(); j++)
    {
        new_E += -Epar.g * polymer[j].r[2] * polymer[j].t[0];
    }

    if (satisfy_self_avoiding_condition_by_group(bead_i, L + 1) && rand_uni(gen) < std::exp(-beta * (new_E - old_E)))
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

int semiflexible_polymer::update_bead_pivot_left(int bead_i)
{
    std::vector<bead> old_beads{}; // starting with bead_i
    double old_E = 0;
    double new_E = 0;
    double old_Tb = 0;
    double new_Tb = 0;
    for (int j = 0; j <= bead_i; j++)
    {
        old_beads.push_back(polymer[j]); // (i),(i+1)... up to (L-1),(L)
    }

    // for any bead_i with a t, which is [0,L-1]
    // everything j < bead_i us updated, since t_i rotate
    // what updates:
    // t(0), t(1), ...  t(i-1)
    // Energ change related
    // bending (i-1)--(i)--(i+1)-- on (i) if i>1 else, no bending changed
    // streching r(i) - r(0)
    // shearing all j<bead_i
    // TODO: wrote up to here continue!!!
    if (bead_i == 0)
    {
        std::cout << "Error: bead_i == 0" << std::endl;
        return 0;
    }
    // note:  bead_i > 0
    old_Tb = calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i].t);
    old_E = 0.5 * Epar.kappa * old_Tb;
    old_E += -Epar.f * (polymer[L].r[0] - polymer[0].r[0]);
    for (int j = 0; j < bead_i; j++)
    {
        old_E += -Epar.g * polymer[j].r[2] * polymer[j].t[0];
    }

    // update
    // 1. find random unit uni_v vector not parallel to t
    std::vector<double> uni_v = rand_uni_vec();
    while (std::abs(inner_product(uni_v, polymer[bead_i - 1].t)) < 1e-6)
    { // uni_v must be not parallel to t
        uni_v = rand_uni_vec();
    }
    // 2. transform uni_v to vector uni_w perpendicular to t
    std::vector<double> uni_w = cross_product(polymer[bead_i - 1].t, uni_v);
    double uni_w_norm = std::sqrt(uni_w[0] * uni_w[0] + uni_w[1] * uni_w[1] + uni_w[2] * uni_w[2]);
    uni_w[0] /= uni_w_norm;
    uni_w[1] /= uni_w_norm;
    uni_w[2] /= uni_w_norm;

    // also need uni_v perpendicular to both w and t
    uni_v = cross_product(polymer[bead_i - 1].t, uni_w);
    double uni_v_norm = std::sqrt(uni_v[0] * uni_v[0] + uni_v[1] * uni_v[1] + uni_v[2] * uni_v[2]);
    uni_v[0] /= uni_v_norm;
    uni_v[1] /= uni_v_norm;
    uni_v[2] /= uni_v_norm;

    // 3. rotate every r(ij) towards uni_w for random angle theta such that cos(theta) uniform in [cos d_theta,1]
    // due to the length of rotated polymer, here use d_theta/(L-bead_i), to adjust acceptance rate
    double d_theta_in_use = (Epar.f == 0 && Epar.g == 0) ? d_theta : d_theta / (1 + (Epar.f + Epar.g) * bead_i);
    double cos_theta = 1 - rand_uni(gen) * (1 - std::cos(1.0 * d_theta_in_use));
    double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
    std::vector<double> r_ij(3, 0);
    double r_ij_norm = 0;
    double r_ij_norm_new = 0;
    double r_ij_v = 0;
    double r_ij_w = 0;
    double r_ij_t = 0;

    double tj_v = 0;
    double tj_w = 0;
    double tj_t = 0;
    double tj_norm = 0;
    // note: txw = v
    // check othorgonal of v,w,t
    if (std::abs(inner_product(uni_v, uni_w)) > 1e-6 || std::abs(inner_product(uni_v, polymer[bead_i - 1].t)) > 1e-6 || std::abs(inner_product(uni_w, polymer[bead_i - 1].t)) > 1e-6)
    {
        std::cout << "Error: uni_v, uni_w, t not orthogonal" << std::endl;
        std::cout << inner_product(uni_v, uni_w) << ":" << inner_product(uni_v, polymer[bead_i].t) << ":" << inner_product(uni_w, polymer[bead_i].t) << std::endl;
        return 0;
    }
    std::vector<double> ti_old = polymer[bead_i].t;
    for (int j = 0; j < bead_i; j++) // polymer[L] has no t
    {
        // rotate tj towards uni_w by theta
        tj_v = inner_product(polymer[j].t, uni_v);
        tj_w = inner_product(polymer[j].t, uni_w);
        tj_t = inner_product(polymer[j].t, ti_old);
        polymer[j].t = {(cos_theta * tj_t - sin_theta * tj_w) * ti_old[0] + (sin_theta * tj_t + cos_theta * tj_w) * uni_w[0] + tj_v * uni_v[0],
                        (cos_theta * tj_t - sin_theta * tj_w) * ti_old[1] + (sin_theta * tj_t + cos_theta * tj_w) * uni_w[1] + tj_v * uni_v[1],
                        (cos_theta * tj_t - sin_theta * tj_w) * ti_old[2] + (sin_theta * tj_t + cos_theta * tj_w) * uni_w[2] + tj_v * uni_v[2]};
        // renormalize t
        tj_norm = std::sqrt(polymer[j].t[0] * polymer[j].t[0] + polymer[j].t[1] * polymer[j].t[1] + polymer[j].t[2] * polymer[j].t[2]);
        polymer[j].t[0] /= tj_norm;
        polymer[j].t[1] /= tj_norm;
        polymer[j].t[2] /= tj_norm;
    }

    // reconstruct r
    for (int k = bead_i - 1; k >= 0; k--)
    {
        polymer[k].r = {polymer[k + 1].r[0] - polymer[k].t[0],
                        polymer[k + 1].r[1] - polymer[k].t[1],
                        polymer[k + 1].r[2] - polymer[k].t[2]};
    }

    // calculate updated energy
    new_Tb = calc_bending_of_two_t(polymer[bead_i - 1].t, polymer[bead_i].t);
    new_E = 0.5 * Epar.kappa * new_Tb;
    new_E += -Epar.f * (polymer[L].r[0] - polymer[0].r[0]);
    for (int j = 0; j < bead_i; j++)
    {
        new_E += -Epar.g * polymer[j].r[2] * polymer[j].t[0];
    }

    if (satisfy_self_avoiding_condition_by_group(0, bead_i) && rand_uni(gen) < std::exp(-beta * (new_E - old_E)))
    {
        E_sys += new_E - old_E;
        Tb_sys += new_Tb - old_Tb;
        return 1;
    }
    else
    {
        for (int j = 0; j <= bead_i; j++)
        {
            polymer[j] = old_beads[j];
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

int semiflexible_polymer::satisfy_self_avoiding_condition_by_group(int bead_i, int bead_j)
{
    // check overlapping between j<bead_i and j>bead_i
    if (bead_i > bead_j)
    {
        // swap, make sure i is before j
        int bead = bead_i;
        bead_i = bead_j;
        bead_j = bead;
    }
    if (bead_j - bead_i < 2)
    {
        // nothing in between
        return 1;
    }
    double rjk_norm2 = 0;

    if (bead_j == L + 1)
    {
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
    }
    else
    {
        for (int j = bead_i + 1; j < bead_j; j++)
        {
            for (int k = 0; k < polymer.size(); k++)
            {
                if (k > bead_i && k < bead_j)
                {
                    // within the group, no need to check
                    continue;
                }
                rjk_norm2 = 0;
                for (int l = 0; l < 3; l++)
                {
                    rjk_norm2 += (polymer[j].r[l] - polymer[k].r[l]) * (polymer[j].r[l] - polymer[k].r[l]);
                }
                if (rjk_norm2 < 1)
                {
                    return 0;
                }
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
            f << "E(energy),Tb(total bending), X(end to end X distance), Y(end to end Y distance), Z(end to end Z distance), Xsign(Z), Zsign(X), R(end to end distance), R^2,"
              << "Rg^2(radius of gyration square), Sxx(gyration tensor component), Syy, Szz, Sxy, Sxz, Syz ";
            for (int i = 0; i < number_of_polymer; i++)
            {
                f << "\n"
                  << obs_ensemble[i].E << "," << obs_ensemble[i].Tb
                  << "," << obs_ensemble[i].X << "," << obs_ensemble[i].Y
                  << "," << obs_ensemble[i].Z << "," << obs_ensemble[i].XsignZ << "," << obs_ensemble[i].ZsignX
                  << "," << obs_ensemble[i].R << "," << obs_ensemble[i].R2 << "," << obs_ensemble[i].Rg2
                  << "," << obs_ensemble[i].Sxx << "," << obs_ensemble[i].Syy << "," << obs_ensemble[i].Szz
                  << "," << obs_ensemble[i].Sxy << "," << obs_ensemble[i].Sxz << "," << obs_ensemble[i].Syz;
            }
        }
        else
        {
            // only save stats
            // find stats
            // calculate average and standard deviation of E
            double avg_E = 0.0, std_E = 0.0, avg_Tb = 0.0, std_Tb = 0.0;
            double avg_X = 0.0, std_X = 0.0, avg_Y = 0.0, std_Y = 0.0, avg_Z = 0.0, std_Z = 0.0;
            double avg_XsignZ = 0.0, std_XsignZ = 0.0, avg_ZsignX = 0.0, std_ZsignX = 0.0;
            double avg_R = 0.0, std_R = 0.0, avg_R2 = 0.0, std_R2 = 0.0;
            double avg_Rg2 = 0.0, std_Rg2 = 0.0;
            double avg_Sxx = 0.0, std_Sxx = 0.0, avg_Syy = 0.0, std_Syy = 0.0, avg_Szz = 0.0, std_Szz = 0.0;
            double avg_Sxy = 0.0, std_Sxy = 0.0, avg_Sxz = 0.0, std_Sxz = 0.0, avg_Syz = 0.0, std_Syz = 0.0;

            std::vector<double> avg_Sq(obs_ensemble[0].Sq.size(), 0.0);
            std::vector<double> std_Sq(obs_ensemble[0].Sq.size(), 0.0);
            std::vector<std::vector<double>> avg_Sq2D(obs_ensemble[0].Sq2D.size(), std::vector<double>(obs_ensemble[0].Sq2D[0].size(), 0.0));
            std::vector<double> avg_tts(obs_ensemble[0].tts.size(), 0.0);
            std::vector<double> std_tts(obs_ensemble[0].tts.size(), 0.0);

            double M = obs_ensemble.size();
            double sqrt_M = std::sqrt(M);

            // get the average
            for (int i = 0; i < obs_ensemble.size(); i++)
            {
                avg_E += obs_ensemble[i].E;
                avg_Tb += obs_ensemble[i].Tb;
                avg_X += obs_ensemble[i].X;
                avg_Y += obs_ensemble[i].Y;
                avg_Z += obs_ensemble[i].Z;
                avg_XsignZ += obs_ensemble[i].XsignZ;
                avg_ZsignX += obs_ensemble[i].ZsignX;
                avg_R += obs_ensemble[i].R;
                avg_R2 += obs_ensemble[i].R2;
                avg_Rg2 += obs_ensemble[i].Rg2;
                avg_Sxx += obs_ensemble[i].Sxx;
                avg_Syy += obs_ensemble[i].Syy;
                avg_Szz += obs_ensemble[i].Szz;
                avg_Sxy += obs_ensemble[i].Sxy;
                avg_Sxz += obs_ensemble[i].Sxz;
                avg_Syz += obs_ensemble[i].Syz;
                for (int j = 0; j < obs_ensemble[i].Sq.size(); j++)
                {
                    avg_Sq[j] += obs_ensemble[i].Sq[j];
                    avg_tts[j] += obs_ensemble[i].tts[j];
                }
                for (int kx = 0; kx < obs_ensemble[i].Sq2D.size(); kx++)
                {
                    for (int ky = 0; ky < obs_ensemble[i].Sq2D[kx].size(); ky++)
                    {
                        avg_Sq2D[kx][ky] += obs_ensemble[i].Sq2D[kx][ky];
                    }
                }
            }
            avg_E /= M;
            avg_Tb /= M;
            avg_X /= M;
            avg_Y /= M;
            avg_Z /= M;
            avg_XsignZ /= M;
            avg_ZsignX /= M;
            avg_R /= M;
            avg_R2 /= M;
            avg_Rg2 /= M;
            avg_Sxx /= M;
            avg_Syy /= M;
            avg_Szz /= M;
            avg_Sxy /= M;
            avg_Sxz /= M;
            avg_Syz /= M;

            if (obs_ensemble[0].Sq.size() != obs_ensemble[0].tts.size())
            {
                std::cout << "Error: Sq and tts size not match" << std::endl;
                return;
            }

            for (int j = 0; j < obs_ensemble[0].Sq.size(); j++)
            {
                avg_Sq[j] /= M;
                avg_tts[j] /= M;
            }
            for (int kx = 0; kx < obs_ensemble[0].Sq2D.size(); kx++)
            {
                for (int ky = 0; ky < obs_ensemble[0].Sq2D[kx].size(); ky++)
                {
                    avg_Sq2D[kx][ky] /= M;
                }
            }

            // get std
            for (int i = 0; i < obs_ensemble.size(); i++)
            {
                std_E += (obs_ensemble[i].E - avg_E) * (obs_ensemble[i].E - avg_E);
                std_Tb += (obs_ensemble[i].Tb - avg_Tb) * (obs_ensemble[i].Tb - avg_Tb);
                std_X += (obs_ensemble[i].X - avg_X) * (obs_ensemble[i].X - avg_X);
                std_Y += (obs_ensemble[i].Y - avg_Y) * (obs_ensemble[i].Y - avg_Y);
                std_Z += (obs_ensemble[i].Z - avg_Z) * (obs_ensemble[i].Z - avg_Z);
                std_XsignZ += (obs_ensemble[i].XsignZ - avg_XsignZ) * (obs_ensemble[i].XsignZ - avg_XsignZ);
                std_ZsignX += (obs_ensemble[i].ZsignX - avg_ZsignX) * (obs_ensemble[i].ZsignX - avg_ZsignX);
                std_R += (obs_ensemble[i].R - avg_R) * (obs_ensemble[i].R - avg_R);
                std_R2 += (obs_ensemble[i].R2 - avg_R2) * (obs_ensemble[i].R2 - avg_R2);
                std_Rg2 += (obs_ensemble[i].Rg2 - avg_Rg2) * (obs_ensemble[i].Rg2 - avg_Rg2);
                std_Sxx += (obs_ensemble[i].Sxx - avg_Sxx) * (obs_ensemble[i].Sxx - avg_Sxx);
                std_Syy += (obs_ensemble[i].Syy - avg_Syy) * (obs_ensemble[i].Syy - avg_Syy);
                std_Szz += (obs_ensemble[i].Szz - avg_Szz) * (obs_ensemble[i].Szz - avg_Szz);
                std_Sxy += (obs_ensemble[i].Sxy - avg_Sxy) * (obs_ensemble[i].Sxy - avg_Sxy);
                std_Sxz += (obs_ensemble[i].Sxz - avg_Sxz) * (obs_ensemble[i].Sxz - avg_Sxz);
                std_Syz += (obs_ensemble[i].Syz - avg_Syz) * (obs_ensemble[i].Syz - avg_Syz);

                for (int j = 0; j < obs_ensemble[i].Sq.size(); j++)
                {
                    std_Sq[j] += (obs_ensemble[i].Sq[j] - avg_Sq[j]) * (obs_ensemble[i].Sq[j] - avg_Sq[j]);
                    std_tts[j] += (obs_ensemble[i].tts[j] - avg_tts[j]) * (obs_ensemble[i].tts[j] - avg_tts[j]);
                }
            }
            std_E = std::sqrt(std_E / M);
            std_Tb = std::sqrt(std_Tb / M);
            std_X = std::sqrt(std_X / M);
            std_Y = std::sqrt(std_Y / M);
            std_Z = std::sqrt(std_Z / M);
            std_XsignZ = std::sqrt(std_XsignZ / M);
            std_ZsignX = std::sqrt(std_ZsignX / M);
            std_R = std::sqrt(std_R / M);
            std_R2 = std::sqrt(std_R2 / M);
            std_Rg2 = std::sqrt(std_Rg2 / M);
            std_Sxx = std::sqrt(std_Sxx / M);
            std_Syy = std::sqrt(std_Syy / M);
            std_Szz = std::sqrt(std_Szz / M);
            std_Sxy = std::sqrt(std_Sxy / M);
            std_Sxz = std::sqrt(std_Sxz / M);
            std_Syz = std::sqrt(std_Syz / M);

            for (int j = 0; j < obs_ensemble[0].Sq.size(); j++)
            {
                std_Sq[j] = std::sqrt(std_Sq[j] / M);
                std_tts[j] = std::sqrt(std_tts[j] / M);
            }
            // write parameters and stats to the file
            f << "stats,L,kappa,f,g,E,Tb,X,Y,Z,Xsign(Z),Zsign(X),R,R2,Rg2,Sxx,Syy,Szz,Sxy,Sxz,Syz,Sq\n";
            f << "\nmean," << L << "," << Epar.kappa << "," << Epar.f << "," << Epar.g << "," << avg_E << "," << avg_Tb
              << "," << avg_X << "," << avg_Y << "," << avg_Z << "," << avg_XsignZ << "," << avg_ZsignX
              << "," << avg_R << "," << avg_R2 << "," << avg_Rg2 << "," << avg_Sxx << "," << avg_Syy
              << "," << avg_Szz << "," << avg_Sxy << "," << avg_Sxz << "," << avg_Syz;
            for (int j = 0; j < avg_Sq.size(); j++)
            {
                f << "," << avg_Sq[j];
            }

            f << "\nstd/sqrt(number of polymer),NA,NA,NA,NA," << std_E / sqrt_M << "," << std_Tb / sqrt_M
              << "," << std_X / sqrt_M << "," << std_Y / sqrt_M
              << "," << std_Z / sqrt_M << "," << std_XsignZ / sqrt_M << "," << std_ZsignX / sqrt_M
              << "," << std_R / sqrt_M << "," << std_R2 / sqrt_M << "," << std_Rg2 / sqrt_M
              << "," << std_Sxx / sqrt_M << "," << std_Syy / sqrt_M << "," << std_Szz / sqrt_M
              << "," << std_Sxy / sqrt_M << "," << std_Sxz / sqrt_M << "," << std_Syz / sqrt_M;
            for (int j = 0; j < std_Sq.size(); j++)
            {
                f << "," << std_Sq[j] / sqrt_M;
            }

            f << "\n qB,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < obs_ensemble[0].qB.size(); j++)
            {
                f << "," << obs_ensemble[0].qB[j];
            }

            f << "\n tts,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < avg_tts.size(); j++)
            {
                f << "," << avg_tts[j];
            }
            f << "\nstd_tts/sqrt(number of polymer),NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < std_tts.size(); j++)
            {
                f << "," << std_tts[j] / sqrt_M;
            }
            f << "\n r(for tts),NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int j = 0; j < obs_ensemble[0].spB.size(); j++)
            {
                f << "," << obs_ensemble[0].spB[j];
            }
            for (int kx = 0; kx < obs_ensemble[0].Sq2D.size(); kx++)
            {
                f << "\nSq2D,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
                for (int ky = 0; ky < obs_ensemble[0].Sq2D[kx].size(); ky++)
                {
                    f << "," << avg_Sq2D[kx][ky];
                }
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
    obs.R2 = R[0] * R[0] + R[1] * R[1] + R[2] * R[2];
    obs.X = R[0];
    obs.Y = R[1];
    obs.Z = R[2];
    obs.XsignZ = (R[2] > 0) ? R[0] : -R[0];
    obs.ZsignX = (R[0] > 0) ? R[2] : -R[2];
    obs.Rg2 = calc_radius_of_gyration_square();

    std::vector<double> Sij = {0, 0, 0, 0, 0, 0};
    Sij = calc_gyration_tensor();
    obs.Sxx = Sij[0];
    obs.Syy = Sij[1];
    obs.Szz = Sij[2];
    obs.Sxy = Sij[3];
    obs.Sxz = Sij[4];
    obs.Syz = Sij[5];

    double qB_i = -50.0 / L * M_PI;                // 0.2*M_PI/L; //0.1/L; ;
    double dqB = 100.0 / L * M_PI / (bin_num - 1); // M_PI;//100.0/L; //M_PI;
    obs.qB.resize(bin_num);
    for (int k = 0; k < bin_num; k++)
    {
        obs.qB[k] = qB_i + dqB * k;
    }
    obs.Sq.resize(bin_num);
    obs.Sq2D = calc_structure_factor_2d(obs.qB);
    /*
    for (int k = 0; k < bin_num; k++)
    {
        obs.qB[k] = qB_i * std::pow(qB_f / qB_i, 1.0 * k / (bin_num - 1)); // uniform in log scale
    }
    obs.Sq = calc_structure_factor(obs.qB);
    */

    obs.spB.resize(bin_num);
    obs.tts.resize(bin_num);
    for (int k = 0; k < bin_num; k++)
    {
        obs.spB[k] = k; // 0, 1, 2, ..
    }
    obs.tts = calc_tangent_pair_correlation(obs.spB);

    return obs;
}

std::vector<double> semiflexible_polymer::calc_structure_factor(std::vector<double> qB)
{
    // measure structure factor
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

std::vector<std::vector<double>> semiflexible_polymer::calc_structure_factor_2d(std::vector<double> qB)
{
    // measrue 2d structure factor
    int bin_num = qB.size();
    int bin0 = (bin_num - 1) / 2;
    if (bin_num % 2 == 0)
    {
        std::cout << "Error: bin_num should be odd" << std::endl;
    }
    if (std::abs(qB[bin0]) > 1e-6)
    {
        std::cout << "Error: qB[bin0] should be 0" << std::endl;
        std::cout << "qB[bin0]" << qB[bin0] << std::endl;
    }

    std::vector<std::vector<double>> R_all{}; // all scattering point's R[axis] value, including all scattering points in each segment

    for (int i = 0; i < size(polymer); i++)
    {
        R_all.push_back(polymer[i].r);
    }
    int N = R_all.size(); // total number of scattering points
    std::vector<double> r{0, 0, 0};
    double qx, qz, SqB_Re_buff, SqB_Im_buff;
    // std::vector<std::vector<double>> SqB(bin_num, std::vector<double>(bin_num, 1.0 / N)); // initialize with 1 due to self overlaping term (see S.H. Chen 1986 eq 18)
    std::vector<std::vector<double>> SqB(bin_num, std::vector<double>(bin_num, 0.0));    // initialize with 1 due to self overlaping term (see S.H. Chen 1986 eq 18)
    std::vector<std::vector<double>> SqB_Re(bin_num, std::vector<double>(bin_num, 0.0)); // real part
    N = L + 1;
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            r[0] = polymer[i].r[0] - polymer[j].r[0];
            r[1] = polymer[i].r[1] - polymer[j].r[1];
            r[2] = polymer[i].r[2] - polymer[j].r[2];
            // calculate S_q
            for (int kx = 0; kx < bin_num; kx++)
            {
                qx = qB[kx];
                for (int kz = bin0; kz < bin_num; kz++)
                {
                    qz = qB[kz];
                    SqB_Re_buff = 2.0 / (N * N) * std::cos(qx * r[0] + qz * r[2]);
                    SqB_Re[kx][kz] += SqB_Re_buff;

                    if (kz != bin0)
                    {
                        SqB_Re[2 * bin0 - kx][2 * bin0 - kz] += SqB_Re_buff;
                    }
                }
            }
        }
    }

    // calculate SqB
    for (int kx = 0; kx < bin_num; kx++)
    {
        for (int kz = 0; kz < bin_num; kz++)
        {
            SqB[kx][kz] = (1.0 / N + SqB_Re[kx][kz]); // * SqB_Re[kx][ky] + SqB_Im[kx][ky] * SqB_Im[kx][ky];
        }
    }
    return SqB;
}

std::vector<double> semiflexible_polymer::calc_tangent_pair_correlation(std::vector<double> spB)
{
    int nbin = spB.size();
    std::vector<double> tts(nbin, 0);
    tts[0] = 1.0;                         // self correlation
    std::vector<int> pair_count(nbin, 0); // count number of t-t pair added to a s
    for (int i = 0; i < L - 1; i++)
    {
        for (int j = i + 1; j < L; j++)
        {
            if (j - i >= nbin)
            {
                continue;
            }
            tts[j - i] += inner_product(polymer[i].t, polymer[j].t);
            pair_count[j - i]++;
        }
    }
    for (int i = 1; i < nbin; i++)
    {
        tts[i] /= pair_count[i];
    }
    return tts;
}

std::vector<double> semiflexible_polymer::calc_gyration_tensor()
{
    double Sxx = 0, Syy = 0, Szz = 0, Sxy = 0, Sxz = 0, Syz = 0;
    int N = polymer.size();
    for (int i = 0; i < polymer.size() - 1; i++)
    {
        for (int j = i; j < polymer.size(); j++)
        {
            Sxx += (polymer[i].r[0] - polymer[j].r[0]) * (polymer[i].r[0] - polymer[j].r[0]);
            Syy += (polymer[i].r[1] - polymer[j].r[1]) * (polymer[i].r[1] - polymer[j].r[1]);
            Szz += (polymer[i].r[2] - polymer[j].r[2]) * (polymer[i].r[2] - polymer[j].r[2]);
            Sxy += (polymer[i].r[0] - polymer[j].r[0]) * (polymer[i].r[1] - polymer[j].r[1]);
            Sxz += (polymer[i].r[0] - polymer[j].r[0]) * (polymer[i].r[2] - polymer[j].r[2]);
            Syz += (polymer[i].r[1] - polymer[j].r[1]) * (polymer[i].r[2] - polymer[j].r[2]);
        }
    }
    Sxx /= N * N;
    Syy /= N * N;
    Szz /= N * N;
    Sxy /= N * N;
    Sxz /= N * N;
    Syz /= N * N;
    return {Sxx, Syy, Szz, Sxy, Sxz, Syz};
}

double semiflexible_polymer::calc_radius_of_gyration_square()
{
    // calculate radius of gyration
    double Rg2 = 0;
    // find center of mass
    double x_c = 0, y_c = 0, z_c = 0;
    for (int i = 0; i < polymer.size(); i++)
    {
        x_c = x_c + polymer[i].r[0];
        y_c = y_c + polymer[i].r[1];
        z_c = z_c + polymer[i].r[2];
    }
    x_c = x_c / L;
    y_c = y_c / L;
    z_c = z_c / L;

    // calculate Rg
    for (int i = 0; i < polymer.size(); i++)
    {
        Rg2 += std::pow(polymer[i].r[0] - x_c, 2) + std::pow(polymer[i].r[1] - y_c, 2) + std::pow(polymer[i].r[2] - z_c, 2);
    }
    Rg2 /= polymer.size();
    return Rg2;
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

double semiflexible_polymer::BesselJ0(double x)
{
    // series expansion for J_0
    double fct = 1;
    double sum = 0;
    for (int k = 0; k < 6; fct *= ++k)
    {
        sum += std::pow(-1, k) * std::pow(x / 2, 2 * k) / std::pow(fct, 2);
        // std::cout << "sum = " << sum << '\n';
    }
    return sum;
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

double semiflexible_polymer::run_MC_sweep(int step_per_sweep)
{
    int bead_i, bead_j, bead_ij;
    double acceptance_rate = 0;
    // int mid_point = int(L / 2);
    for (int n = 0; n < step_per_sweep; n++)
    {
        bead_ij = 2 + int(rand_uni(gen) * Epar.kappa); //~[2,L]
        bead_i = int(rand_uni(gen) * (L - bead_ij));
        bead_i += bead_ij;
        bead_j = bead_i + bead_ij;

        acceptance_rate += update_bead_crankshaft(bead_i, bead_j);
        /*
        // take between [0,mid_poin-bead_ij] and [mid_point,L-bead_ij]
        if (bead_i > mid_point - bead_ij)
        {
            bead_i += bead_ij;
            bead_j = bead_i + bead_ij;
            bead_j = std::min(bead_j, L - 1);
        }
        else
        {
            bead_j = bead_i + bead_ij;
            bead_j = std::min(bead_j, mid_point);
        }
        */

        bead_i = int(rand_uni(gen) * L); // take from [1,L-1]
        update_bead_pivot_right(bead_i);
    }
    return acceptance_rate / step_per_sweep;
}

void semiflexible_polymer::run_simultion(int therm_sweep, int MC_sweeps, int step_per_sweep, std::string folder, std::string finfo, int bin_num, int save_more_config)
{
    std::vector<observable> obs_ensemble;

    double acceptance_rate = 0;

    beta = 0; // randomization
    for (int i = 0; i < therm_sweep; i++)
    {
        std::cout << "randomization(beta=0) sweep " << i << " out of " << therm_sweep << " (" << (i * 100) / therm_sweep << "%)\r";
        run_MC_sweep(step_per_sweep);
    }

    std::cout << "\n";
    // thermalization using tempering
    for (int i = 0; i < therm_sweep; i++)
    {
        beta = (i + 1) * 1.0 / therm_sweep; // tempering
        std::cout << "thermalization/tempering (beta -> 1) sweep " << i << " out of " << therm_sweep << " (" << (i * 100) / therm_sweep << "%)\r";
        run_MC_sweep(step_per_sweep);
    }

    // data collection run
    // TODO: consider implementing irreversible algorithem like event-chain here if too slow
    std::cout << "\n";
    if (beta != 1)
    {
        std::cout << "Error: beta is not 1 at the end of thermalization\n";
        beta = 1;
    }
    if (save_more_config)
    {
        std::cout << "saving more config, creating foler: " << folder + "/" + finfo << std::endl;
        std::filesystem::create_directory(folder + "/" + finfo);
    }

    for (int i = 0; i < MC_sweeps; i++)
    {
        std::cout << "MC sweep " << i << " out of " << MC_sweeps << " (" << (i * 100) / MC_sweeps << "%)\r";
        acceptance_rate += run_MC_sweep(step_per_sweep);

        obs_ensemble.push_back(measure_observable(bin_num));
        if (i % 100 == 0 && save_more_config)
        {
            save_polymer_to_file(folder + "/" + finfo + "/config_" + std::to_string(int(i / 100)) + ".csv");
        }
    }
    std::cout << "\n";
    std::cout << "crankshaft_acceptance rate:" << acceptance_rate / MC_sweeps << std::endl;

    save_polymer_to_file(folder + "/config_" + finfo + ".csv");
    // save_observable_to_file(folder + "/obs_MC_" + finfo + ".csv", obs_ensemble, true);
    save_observable_to_file(folder + "/obs_" + finfo + ".csv", obs_ensemble, false);
}