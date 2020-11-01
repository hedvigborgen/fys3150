#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H

class VelocityVerlet
{
public:
    double m_dt;
    VelocityVerlet(double dt);
    void integrateOneStep(class SolarSystem &system, double beta);
};

#endif // VELOCITYVERLET_H
