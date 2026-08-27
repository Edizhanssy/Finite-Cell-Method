//
// Created by Edizhan Yildiz on 06.08.26.
//

#ifndef FINITE_CELL_METHOD_CONFIG_H
#define FINITE_CELL_METHOD_CONFIG_H

#pragma once
#include <array>
#include <vector>

namespace fcm {

    struct DomainSpan {
        double x0;
        double x1;
        int material;          // 0 = physical, 1 = fictitious
    };

    struct Config {
        double L          = 3.0;
        int    n_elements = 2;
        double E          = 1.0;
        double A          = 1.0;
        double alpha      = 1e-8;
        int    p          = 15;
        int    max_depth  = 10;
        double penalty    = 1e5;
        double disp_load  = -1.0;

        std::vector<DomainSpan> domains{
            {0.0, 1.0, 0}, {1.0, 7.0 / 3.0, 1}, {7.0 / 3.0, 3.0, 0}};

        double load_span0 = 0.0;
        double load_span1 = 1.0;
        double load_amp   = 1.0 / 20.0;
        double load_freq  = 4.0 * 3.14159265358979323846;

        int n_gauss() const { return p + 1; }
        int n_modes() const { return p + 1; }

        int    domain_index(double x) const;
        double material_factor(int material) const;
        double body_load(double x) const;
        double to_local(double x) const { return (2.0 * x - L) / L; }
        std::array<double, 2> fictitious_span() const;
    };

}


#endif //FINITE_CELL_METHOD_CONFIG_H