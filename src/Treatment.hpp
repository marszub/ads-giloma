#include <cmath>

class Treatment
{
public:
    Treatment() = default;
    Treatment(double absorptionRate, double decayRate, double dose,
              double firstDoseTime, int dosesNum,
              double timeBetweenDoses)
        : absorptionRate(absorptionRate), decayRate(decayRate),
          dose(dose), firstDoseTime(firstDoseTime), dosesNum(dosesNum),
          timeBetweenDoses(timeBetweenDoses) {}

    double operator()(double x, double y, double t) const
    {
        double f = 0;
        for (int i = 0; i < dosesNum; i++)
        {
            double peak_center = firstDoseTime + i * timeBetweenDoses;
            double absorptionPhase =
                dose * exp(-absorptionRate * pow((t - peak_center), 2)) *
                (t <= peak_center);
            double decayPhase =
                dose * exp(-decayRate * pow((t - peak_center), 2)) *
                (t > peak_center);
            f += absorptionPhase + decayPhase;
        }
        return f;
    }

private:
    double absorptionRate{0};
    double decayRate{0};
    double dose{0};
    double firstDoseTime{0};
    int dosesNum{0};
    double timeBetweenDoses{0};
};