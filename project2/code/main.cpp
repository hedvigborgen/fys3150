#include "project2.hpp"

int main(int argc, char const *argv[]) {
  methods object;
  object.n = 3; //atoi(argv[1]);
  object.d = 1; //atoi(argv[2]);
  object.a = 2; //atoi(argv[3]);

  object.A = mat(object.n,object.n);

  for (int i = 1; i < object.n; i++){
    object.A(i,i) = object.d;
    object.A(i-1,i) = object.a;
    object.A(i,i-1) = object.a;
   }

   object.diagonalize(object.eigval, object.eigvec, object.A);



  /*
  object.B = mat(n,n).eye;


  double eps = pow(10,-8);
  object.offA = 100;

  while (offA > eps){
    object.jacobi(A)
    object.offA(A)
  }
*/
  return 0;
}
