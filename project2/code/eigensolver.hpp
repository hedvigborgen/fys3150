#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Eigensolver{
  private:
    double m_max_off_d;
    vec m_init_eigval;
    mat m_init_eigvec;
    mat m_eigvec;
    double m_V;
    int m_count;

  public:
    int m_n;
    mat m_A;
    vec m_eigval;
    int m_k;
    int m_l;
    //void V(double rho, double omega);
    void initialize(double rho_max);
    void rotation(mat A);
    void max_val(mat A);
    void diagonalize(mat A);
    void print();
    void print_test();

};