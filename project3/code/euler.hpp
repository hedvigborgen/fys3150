#ifndef EULER_H
#define EULER_H

class Euler
{
public:
    double m_dt;
    Euler(double dt);
    void integrateOneStep(class SolarSystem &system, double beta);
};

#endif // EULER_H
