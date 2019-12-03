#include <vector>
#include <cmath>

class LevelSymmetric
{
public:
  LevelSymmetric(const unsigned int& N, const double& mu_1) : _N(N), _mu_1(mu_1) {
    calculateMuVector();
    calculateQuadratureAndWeight();
  }
  const std::vector<double> & getWeight() const {
      return _weights;
  }
  const std::vector<std::vector<double>> & getQuadrature() const {
      return _quad_points;
  }
  const std::vector<double> & getMu() const {
      return _mu;
  }

private:
  unsigned int _N;
  double _mu_1;
  std::vector<double> _weights, _mu;
  std::vector<std::vector<double>> _quad_points;
  void calculateMuVector();
  void calculateQuadratureAndWeight();
};

void LevelSymmetric::calculateMuVector(){
  _mu.push_back(_mu_1);
  for (auto i = 2; i <= _N/2; i++) {
    _mu.push_back(sqrt(_mu_1*_mu_1 + (i-1)*2*(1-3*_mu_1*_mu_1)/(_N-2)));
  }
}

void LevelSymmetric::calculateQuadratureAndWeight(){
  for (auto i = 1; i <= _N/2; i++) {
    for (auto j = 1; j <= _N/2-i+1; j++) {
      _quad_points.push_back( {_mu[i-1], _mu[j-1], _mu[_N/2+2-i-j-1]} );
    }
  }
  std::vector<unsigned int> temp_vec;
  if(_N == 8) temp_vec={1,2,2,1,2,3,2,2,2,1};
  for (auto i = 0; i < temp_vec.size(); i++) {
    if(temp_vec[i]==1) _weights.push_back(0.1209877);
    else if(temp_vec[i]==2) _weights.push_back(0.0907407);
    else if(temp_vec[i]==3) _weights.push_back(0.0925926);
  }
}
