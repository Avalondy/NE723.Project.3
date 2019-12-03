/*
This c++ code for LU decomposition was obtained from obtained from wikipedia with some minor changes
https://en.wikipedia.org/wiki/LU_decomposition
*/
#include <vector>
#include <cmath>

/* INPUT: _A - array of pointers to rows of _A square matrix having dimension _N
 *        _tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix _A is changed, it contains _A copy of both matrices L-E and U as _A=(L-E)+U such that _P*_A=L*U.
 *        The permutation matrix is not stored as _A matrix, but in an integer vector _P of size _N+1
 *        containing column indexes where the permutation matrix has "1". The last element _P[_N]=S+_N,
 *        where S is the number of row exchanges needed for determinant computation, det(_P)=(-1)^S
 */
class LUdecomposition
{
public:
  LUdecomposition(const std::vector<std::vector<double>>& A, const std::vector<double>& b) : _A(A), _b(b){
    _N=_A.size();
    _P.resize(_N+1);
  }
  std::vector<double> LUPSolve();
private:
  std::vector<std::vector<double>> _A;
  std::vector<double> _b;
  std::vector<int> _P;
  int _N;
  double _tol=1e-8;

  int LUPDecompose();
};

int LUdecomposition::LUPDecompose()
{
  int j, imax;
  double maxA, absA;
  std::vector<double> ptr;

  for (auto i = 0; i <= _N; i++)
    _P[i] = i; //Unit permutation matrix, _P[_N] initialized with _N
  for (auto i = 0; i < _N; i++) {
    maxA = 0.0;
    imax = i;
    for (auto k = i; k < _N; k++)
      if ((absA = std::abs(_A[k][i])) > maxA) {
          maxA = absA;
          imax = k;
      }
    if (maxA < _tol) return 0; //failure, matrix is degenerate
    if (imax != i) {
        //pivoting _P
        j = _P[i];
        _P[i] = _P[imax];
        _P[imax] = j;
        //pivoting rows of _A
        ptr = _A[i];
        _A[i] = _A[imax];
        _A[imax] = ptr;

        //counting pivots starting from _N (for determinant)
        _P[_N]++;
    }
    for (j = i + 1; j < _N; j++) {
      _A[j][i] /= _A[i][i];
      for (auto k = i + 1; k < _N; k++)
          _A[j][k] -= _A[j][i] * _A[i][k];
    }
  }
  return 1;  //decomposition done
}

/* INPUT: _A,_P filled in LUPDecompose; b - rhs vector; _N - dimension
 * OUTPUT: x - solution vector of _A*x=b
 */
std::vector<double> LUdecomposition::LUPSolve() {
  LUPDecompose();
  std::vector<double> x(_N);
  for (auto i = 0; i < _N; i++) {
      x[i] = _b[_P[i]];
      for (auto k = 0; k < i; k++)
          x[i] -= _A[i][k] * x[k];
  }
  for (auto i = _N - 1; i >= 0; i--) {
      for (auto k = i + 1; k < _N; k++)
          x[i] -= _A[i][k] * x[k];
      x[i] = x[i] / _A[i][i];
  }
  return x;
}

/* INPUT: _A,_P filled in LUPDecompose; _N - dimension
 * OUTPUT: _IA is the inverse of the initial matrix
 */
// void LUPInvert(double **_A, int *_P, int _N, double **_IA)
// {
//   for (int j = 0; j < _N; j++) {
//     for (int i = 0; i < _N; i++) {
//         if (_P[i] == j) _IA[i][j] = 1.0;
//         else            _IA[i][j] = 0.0;
//         for (int k = 0; k < i; k++)
//             _IA[i][j] -= _A[i][k] * _IA[k][j];
//     }
//     for (int i = _N - 1; i >= 0; i--) {
//         for (int k = i + 1; k < _N; k++)
//             _IA[i][j] -= _A[i][k] * _IA[k][j];
//
//         _IA[i][j] = _IA[i][j] / _A[i][i];
//     }
//   }
// }

/* INPUT: _A,_P filled in LUPDecompose; _N - dimension.
 * OUTPUT: Function returns the determinant of the initial matrix
 */
// double LUPDeterminant(double **_A, int *_P, int _N)
// {
//   double det = _A[0][0];
//   for (int i = 1; i < _N; i++)
//       det *= _A[i][i];
//
//   if ((_P[_N] - _N) % 2 == 0)
//       return det;
//   else
//       return -det;
// }
