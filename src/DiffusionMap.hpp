#ifndef DIFFUSION_MAP_HPP
#define DIFFUSION_MAP_HPP

#include <vector>
#include <cmath>
#include "ads/simulation/config.hpp"

class DiffusionMap {
private:
    ads::config_2d config;

public:
    DiffusionMap(ads::config_2d config) : config(config) {}

    double operator()(double x, double y) const {
        double res = 0.0;
        double xDomainBottom = config.x.a;
        double xDomainTop = config.x.b;
        double yDomainBottom = config.y.a;
        double yDomainTop = config.y.b;

        double xDist = (x - xDomainBottom) / (xDomainTop - xDomainBottom) - 0.5;
        double yDist = (y - yDomainBottom) / (yDomainTop - yDomainBottom) - 0.5;
        double distSquared = xDist * xDist + yDist * yDist;
        if (distSquared < 0.49 * 0.49) {
            res = 0.13;
        }
        if (distSquared < 0.15 * 0.15) {
            res = 0.013;
        }
        return res;
    }
};

#endif // DIFFUSION_MAP_HPP
