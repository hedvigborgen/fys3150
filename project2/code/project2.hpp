#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class methods{
  private:
  protected:
  public:
    int n;
    double a;
    double d;
    vec eigval;
    mat eigvec;
    mat A;

    /*
    mat B;
    int k;
    int l;
    double offA;
    vec lambda;


    void jacobi();
    void offA();
    void print();
    */
    void diagonalize();

};
