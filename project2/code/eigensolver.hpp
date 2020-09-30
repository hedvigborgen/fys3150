#include <iostream>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

class Eigensolver{
  private:
    double m_max_off_d;
    vec m_init_eigval;
    mat m_init_eigvec;
    mat m_eigvec;
    vec m_eigval;
    double m_V;
    int m_count;
    vec m_lambda;
    vec m_u;
    mat m_vec;


  public:
    int m_n;
    mat m_A;
    vec m_val;
    int m_k;
    int m_l;
    void initialize(double rho_max, double omega_r);
    void rotation(mat A);
    void max_val(mat A);
    void diagonalize(mat A);
    void print_count();
    void print_eigvals(double omega_r);
    void compare_eigvecs(double rho_max);
    void difference();
    void eigenvecs(double rho_max, double omega_r);

};
