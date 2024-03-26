// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "tumor_2d.hpp"
// #include <torch/torch.h>
// #include <iostream>

int main()
{
    // torch::Tensor tensor = torch::rand({2, 3});
    // std::cout << tensor << std::endl;
    int p = 2;
    int nPoints = 5; // 40
    double spacialDomain = 100;
    int stepsPerDay = 2; // 200
    int days = 100;

    ads::dim_config dim{p, nPoints, 0, spacialDomain};
    ads::timesteps_config steps{days * stepsPerDay, 1.0 / static_cast<double>(stepsPerDay)};
    int derivatives = 1;

    ads::config_2d c{dim, dim, steps, derivatives};
    InitialCondition ic(60.0, 60.0, 0.4, 10);
    DiffusionMap diffusion_map(c);
    Treatment treatment0;
    Treatment treatment2(2.0, 0.02, 0.02, 30.0, 5, 15.0);
    Treatment treatment3(2.0, 0.02, 0.04, 30.0, 5, 15.0);
    ads::problems::tumor_2d sim{c, ic, diffusion_map, treatment3};
    sim.run();
}
