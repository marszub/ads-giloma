#include <cmath>

class InitialCondition {
private:
    double xOrigin, yOrigin;
    double a, b, c;

public:
    InitialCondition(double x, double y, double h, double r) {
        xOrigin = x;
        yOrigin = y;
        a = -h / (2 * pow(r, 2));
        b = -h / (2 * r);
        c = h;
    }

    double operator()(double x, double y) const {
        double d = sqrt(pow(x - xOrigin, 2) + pow(y - yOrigin, 2));
        double res = a * pow(d, 2) + b * d + c;
        res = res > 0 ? res : 0;
        return res;
    }
};
