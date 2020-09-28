#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Eigensolver{
  private:
    vec m_init_eigval;
    mat m_init_eigvec;
    mat m_eigvec;
    double m_V;

  public:
    vec m_eigval;
    mat m_A;
    int m_k;
    int m_l;
    double m_max_off_d;
    void V(double rho, double omega);
    void initialize(int n);
    void rotation(int n, mat A);
    void max_val(int n, mat A);
    void diagonalize(int n, mat A);

};
