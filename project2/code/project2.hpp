#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Eigensolver{
  private:
    vec m_init_eigval;
    mat m_init_eigvec;
    vec m_eigval;
    mat m_eigvec;
    mat m_A;
    double m_V;
    double m_max_off_d;
    int m_k;
    int m_l;

  public:
    void V(double rho, double omega);
    void initialize(int n);
    void rotation(int n);
    void max_len(int n);
    void diagonalize(int n);

};
