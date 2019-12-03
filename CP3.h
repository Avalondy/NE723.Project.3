//header file of CP3.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "RayTracing.h"

enum class TestClass {TestA, TestB};

class TransportEquation
{
public:
  TransportEquation(const TestClass& test_indicator) : _test_indicator(test_indicator){}
  void readInputs(); //as the names suggest
  void initializeParameters();
  void multiGroupSolve();
  void printResults();

private:
  TestClass _test_indicator; //indicate which problem
  SegmentsData _SD_obj;
  PinCellMaterial _PCM_obj;
  int _num_of_lattice_x, _num_of_lattice_y, _num_of_lattice_total, _num_of_cell_per_pin;
  int _N, _G, _M_PO, _M_AZ; //as the names suggest
  double _epsilon_PI, _epsilon_k, _epsilon_transp, _k_init;
  std::vector<double> _q, _k, _radius, _cell_area;
  double _lattice_width, _lattice_height;
  std::vector<double> _rho_PI, _rho_k, _phi_diff_norm, _k_diff_norm;
  std::vector<std::vector<std::vector<double>>> _phi_cell, _psi_bound, _length_modified;

  std::vector<double> sourceIteration(const std::vector<double>& Q, const int& gg);
  void solveTransport(const std::vector<double>&, const int&, const std::vector<double>&);
  std::vector<double> solveSingleLine(const int& m_po, const int& m_az, const int& n, const std::vector<double>& Q, const int& gg, const std::vector<double>& phi_cell);
  void printVector(const std::vector<std::vector<double>>&) const;
  void printVector(const std::vector<double>&) const;
  void printVector(const std::vector<int>&) const;
  std::vector<double> vectorSubtract(std::vector<double>, const std::vector<double>&) const;
  std::vector<std::vector<double>> vectorSubtract(std::vector<std::vector<double>> v1, const std::vector<std::vector<double>>& v2) const;
  double normInf(std::vector<double>) const;
  double normInf(std::vector<std::vector<double>> vec) const;
  void scalarFluxNormalize(std::vector<std::vector<double>>& vec);
};

void TransportEquation::readInputs()
{
  if(_test_indicator == TestClass::TestA) {_num_of_lattice_x=1; _num_of_lattice_y=1;}
  else if(_test_indicator == TestClass::TestB) {_num_of_lattice_x=3; _num_of_lattice_y=3;}
  else throw std::invalid_argument( "Invalid value of test_indicator!" );
  RayTracing RT_obj(_num_of_lattice_x, _num_of_lattice_y);
  _SD_obj = RT_obj.trackSegmentation();
  _PCM_obj = PinCellMaterial(_num_of_lattice_x, _num_of_lattice_y);

  _radius = {0.135,0.27,0.405,0.54,0.585};
  _num_of_lattice_total = _num_of_lattice_x * _num_of_lattice_y;
  _num_of_cell_per_pin = _radius.size()+1;
  _N = _num_of_cell_per_pin*_num_of_lattice_total;
  _M_PO = _SD_obj._polar_theta.size();
  _M_AZ = _SD_obj._azim_phi.size();
  _lattice_width = 1.26;
  _lattice_height = 1.26;
  _G = 7;
  _epsilon_k=1e-5; _epsilon_PI=1e-5; _epsilon_transp=1e-6;
}

void TransportEquation::initializeParameters()
{
  _k_init = 1;
  //cell area
  _cell_area.resize(_num_of_cell_per_pin);
  for (auto p = 0; p < _num_of_cell_per_pin; p++) {
    if(p==0) _cell_area[p] = M_PI*_radius[p]*_radius[p];
    else if(p==_num_of_cell_per_pin-1) _cell_area[p] = _lattice_width*_lattice_height - M_PI*_radius[p-1]*_radius[p-1];
    else _cell_area[p] = M_PI*_radius[p]*_radius[p] - M_PI*_radius[p-1]*_radius[p-1];
  }
  //modified length
  _length_modified = _SD_obj._length;
  for (auto p = 0; p < _N; p++) {
    //calculate f_p
    double temp_denominator = 0;
    for (auto m_po = 0; m_po < _M_PO; m_po++) {
      for (auto m_az = 0; m_az < _M_AZ; m_az++) {
        for (auto n = 0; n < _SD_obj._length[m_az].size(); n++) {
          for (auto l = 0; l < _SD_obj._length[m_az][n].size(); l++) {
            if(_SD_obj._cell[m_az][n][l]==p)
              temp_denominator += _SD_obj._weight_po[m_po]*_SD_obj._weight_az[m_az] * _SD_obj._length[m_az][n][l] * _SD_obj._width[m_az][n];
          }
        }
      }
    }
    temp_denominator = temp_denominator/(4*M_PI);
    double f_p = _cell_area[p%_num_of_cell_per_pin] / temp_denominator;
    //calculate _length_modified
    for (auto m_po = 0; m_po < _M_PO; m_po++) {
      for (auto m_az = 0; m_az < _M_AZ; m_az++) {
        for (auto n = 0; n < _SD_obj._length[m_az].size(); n++) {
          for (auto l = 0; l < _SD_obj._length[m_az][n].size(); l++) {
            _length_modified[m_az][n][l] = f_p * _SD_obj._length[m_az][n][l];
          }
        }
      }
    }
  }
}

void TransportEquation::multiGroupSolve()
{
  _k.push_back(_k_init);
  _phi_cell.push_back(std::vector<std::vector<double>>(_N, std::vector<double>(_G, 1.0) ));
  double temp_while;
  do {
    std::cout << "\nouter iteration = " << _k.size() << '\n' << '\n';
    std::cout << "k" << '\n';
    printVector(_k);
    //calculate q_up and q_fission
    std::vector< std::vector<double> > q_up(_N, std::vector<double>(_G, 0.0)), q_down(_N, std::vector<double>(_G, 0.0)), q_fission(_N, std::vector<double>(_G, 0.0));
    for (auto n = 0; n < _N; n++) {
      for (auto g = 0; g < _G; g++) {
        for (auto gp = g+1; gp < _G; gp++) {
          q_up[n][g] += _PCM_obj._Sigma_s[_SD_obj._cell_to_material[n]][gp][g] * _phi_cell[_k.size()-1][n][gp];
        }
        for (auto gp = 0; gp < _G; gp++) {
          q_fission[n][g] += _PCM_obj._nu_f[_SD_obj._cell_to_material[n]][gp] * _PCM_obj._Sigma_f[_SD_obj._cell_to_material[n]][gp] * _phi_cell[_k.size()-1][n][gp] *_PCM_obj._chi[_SD_obj._cell_to_material[n]][g]/_k.back();
        }
      }
    }

    _phi_cell.push_back(std::vector<std::vector<double>>(_N, std::vector<double>(_G, 0.0)));
    for (auto g = 0; g < _G; g++) {
      // calculate q_down
      for (auto n = 0; n < _N; n++) {
        for (auto gp = 0; gp <= g-1; gp++) {
          q_down[n][g] += _PCM_obj._Sigma_s[_SD_obj._cell_to_material[n]][gp][g] * _phi_cell[_k.size()][n][gp];
        }
      }
      // calculate source to be used in sourceIteration
      std::vector<double> q_g_total(_N, 0.);
      for (auto n = 0; n < _N; n++) {
        q_g_total[n] = q_down[n][g] + q_up[n][g] + q_fission[n][g];
      }
      // sourceIteration solve
      std::vector<double> output_phi_cell = sourceIteration(q_g_total, g);;
      for (auto n = 0; n < _N; n++) _phi_cell[_k.size()][n][g] = output_phi_cell[n];
    }
    // calculate k_eff
    double temp_sum_up=0., temp_sum_down=0.;
    for (auto g = 0; g < _G; g++) {
      for(auto i=0; i<_N; ++i){
        temp_sum_up +=   _cell_area[i%_num_of_cell_per_pin]*_PCM_obj._nu_f[_SD_obj._cell_to_material[i]][g]*_PCM_obj._Sigma_f[_SD_obj._cell_to_material[i]][g]*_phi_cell[_phi_cell.size()-1][i][g];
        temp_sum_down += _cell_area[i%_num_of_cell_per_pin]*_PCM_obj._nu_f[_SD_obj._cell_to_material[i]][g]*_PCM_obj._Sigma_f[_SD_obj._cell_to_material[i]][g]*_phi_cell[_phi_cell.size()-2][i][g];
      }
    }
    _k.push_back(_k.back() * temp_sum_up / temp_sum_down);
    scalarFluxNormalize(_phi_cell.back());
    double temp_norm_up, temp_norm_down;
    if(_k.size()<=2) continue;
    else{
      temp_norm_up   = normInf(vectorSubtract(_phi_cell[_phi_cell.size()-1], _phi_cell[_phi_cell.size()-2]));
      temp_norm_down = normInf(vectorSubtract(_phi_cell[_phi_cell.size()-2], _phi_cell[_phi_cell.size()-3]));
      _rho_PI.push_back(temp_norm_up/temp_norm_down);
      _rho_k.push_back( std::abs(_k[_k.size()-1] - _k[_k.size()-2]) / std::abs(_k[_k.size()-2] - _k[_k.size()-3]) );
      _phi_diff_norm.push_back(temp_norm_up);
      _k_diff_norm.push_back(std::abs(_k[_k.size()-1] - _k[_k.size()-2]));
      temp_while = temp_norm_up;
    }
  } while( (_k.size()<=2) || (temp_while>_epsilon_PI*(1-_rho_PI.back())) || (std::abs(_k[_k.size()-1]-_k[_k.size()-2])>_epsilon_k*(1-_rho_k.back())));
}

std::vector<double> TransportEquation::sourceIteration(const std::vector<double>& Q, const int& gg)
{
  double rho_transp=1.;
  std::vector<std::vector<double>> phi_cell;
  std::vector<double> temp_vec;
  for (auto p = 0; p < _N; p++) temp_vec.push_back(_phi_cell[_k.size()-1][p][gg]);
  phi_cell.push_back(temp_vec);
  do {
    std::cout << "inner iteration = " << phi_cell.size() << '\n';
    // printVector(phi_cell.back());
    solveTransport(Q, gg, phi_cell.back());
    std::vector<double> phi_cell_bar(_N, 0.);
    std::vector<std::vector<std::vector<double>>> psi_cell;
    psi_cell.resize(_M_AZ*_M_PO);
    for (auto m_po = 0; m_po < _M_PO; m_po++) {
      for (auto m_az = 0; m_az < _M_AZ; m_az++) {
        psi_cell[m_po*_M_AZ+m_az].resize(_SD_obj._length[m_az].size());
        for (auto n = 0; n < _SD_obj._length[m_az].size(); n++) {
          psi_cell[m_po*_M_AZ+m_az][n].resize(_SD_obj._length[m_az][n].size(), 0.);
          for (auto l = 0; l < _SD_obj._length[m_az][n].size(); l++) {
            // calculate tau and alpha
            double tau = _PCM_obj._Sigma_t[ _SD_obj._material[m_az][n][l] ][gg]*_length_modified[m_az][n][l]/sin(_SD_obj._polar_theta[m_po]);
            double alpha;
            if(tau<1e-5) alpha = 0.5 - tau/12.;
            else alpha = 1./tau - std::exp(-tau)/(1-std::exp(-tau));
            psi_cell[m_po*_M_AZ+m_az][n][l] = alpha*_psi_bound[m_po*_M_AZ+m_az][n][l] + (1-alpha)*_psi_bound[m_po*_M_AZ+m_az][n][l+1];
          }
        }
      }
    }

    for (auto p = 0; p < _N; p++) {
      double temp_sum = 0;
      for (auto m_po = 0; m_po < _M_PO; m_po++) {
        for (auto m_az = 0; m_az < _M_AZ; m_az++) {
          for (auto n = 0; n < psi_cell[m_az].size(); n++) {
            for (auto l = 0; l < psi_cell[m_az][n].size(); l++) {
              if(_SD_obj._cell[m_az][n][l]==p)
                temp_sum += _SD_obj._weight_po[m_po]*_SD_obj._weight_az[m_az] * _length_modified[m_az][n][l] * _SD_obj._width[m_az][n] * psi_cell[m_po*_M_AZ+m_az][n][l];
            }
          }
        }
      }
      phi_cell_bar[p] = temp_sum/_cell_area[p%_num_of_cell_per_pin];
    }
    phi_cell.push_back(phi_cell_bar);
    if(phi_cell.size() >= 3){
      rho_transp = normInf(vectorSubtract(phi_cell[phi_cell.size()-1], phi_cell[phi_cell.size()-2])) /
                   normInf(vectorSubtract(phi_cell[phi_cell.size()-2], phi_cell[phi_cell.size()-3]));
    }
  } while( normInf(vectorSubtract(phi_cell[phi_cell.size()-1], phi_cell[phi_cell.size()-2])) > _epsilon_transp * (1./rho_transp-1) );
  return (phi_cell.back());
}

void TransportEquation::solveTransport(const std::vector<double>& Q, const int& gg, const std::vector<double>& phi_cell)
{
  //initialize _psi_bound which should have first dimension = M_PO*M_AZ
  _psi_bound.resize(_M_AZ*_M_PO);
  for (auto m_po = 0; m_po < _M_PO; m_po++) {
    for (auto m_az = 0; m_az < _M_AZ; m_az++) {
      _psi_bound[m_po*_M_AZ+m_az].resize(_SD_obj._point[m_az].size());
      for (auto n = 0; n < _SD_obj._point[m_az].size(); n++) {
        _psi_bound[m_po*_M_AZ+m_az][n].resize(_SD_obj._point[m_az][n].size(), 0.);
      }
    }
  }
  // int iter=0;
  //until reflective B.C.s are satisfied.
  for (auto m_po = 0; m_po < _M_PO; m_po++) {
    for (auto m_az = 0; m_az < _M_AZ; m_az++) {
      for (auto n = 0; n < _SD_obj._point[m_az].size(); n++) {
        std::vector<double> start_point = _SD_obj._point[m_az][n][0], end_point;
        int temp_m_az = m_az, temp_n = n, temp_m_az_last = -1, temp_n_last = -1;
        do {
          if(temp_m_az_last != -1 && temp_n_last != -1){
            // apply reflective B.C.
            _psi_bound[m_po*_M_AZ+temp_m_az][temp_n][0] = _psi_bound[m_po*_M_AZ+temp_m_az_last][temp_n_last].back();
          }
          end_point = solveSingleLine(m_po, temp_m_az, temp_n, Q, gg, phi_cell);
          temp_m_az_last = temp_m_az; temp_n_last = temp_n;
          temp_m_az = _SD_obj._reflect_point_m_az[temp_m_az][temp_n];
          temp_n = _SD_obj._reflect_point_n[temp_m_az][temp_n];
        } while( (std::abs(start_point[0]-end_point[0])+std::abs(start_point[1]-end_point[1])) > 1E-12 || std::abs(_psi_bound[m_po*_M_AZ+m_az][n][0]-_psi_bound[m_po*_M_AZ+temp_m_az_last][temp_n_last].back()) > 1E-8);
      }
    }
  }
}

std::vector<double> TransportEquation::solveSingleLine(const int& m_po, const int& m_az, const int& n, const std::vector<double>& Q, const int& gg, const std::vector<double>& phi_cell)
{
  for (auto l = 0; l < _SD_obj._point[m_az][n].size()-1; l++) {
    double exp_term = std::exp( -_PCM_obj._Sigma_t[ _SD_obj._material[m_az][n][l] ][gg]*_length_modified[m_az][n][l]/sin(_SD_obj._polar_theta[m_po]));
    double Q_over_Sigma_t_term;
    if(std::abs(exp_term-1)<1e-12) Q_over_Sigma_t_term = 0.;
    else Q_over_Sigma_t_term = 1/(4*M_PI)*( _PCM_obj._Sigma_s[ _SD_obj._material[m_az][n][l] ][gg][gg] * phi_cell[_SD_obj._cell[m_az][n][l]] + Q[_SD_obj._cell[m_az][n][l]] )
                               /_PCM_obj._Sigma_t[ _SD_obj._material[m_az][n][l] ][gg];
    _psi_bound[m_po*_M_AZ+m_az][n][l+1] = _psi_bound[m_po*_M_AZ+m_az][n][l]*exp_term + Q_over_Sigma_t_term*(1-exp_term);
  }
  return _SD_obj._point[m_az][n].back();
}

void TransportEquation::printResults(){
  std::cout << "iteration #: " << _k.size() << '\n';
  std::cout << "k" << '\n';
  printVector(_k);
  std::cout << "norm of phi difference versus s" << '\n';
  printVector(_phi_diff_norm);
  std::cout << "norm of k difference versus s" << '\n';
  printVector(_k_diff_norm);
  std::cout << "rho_PI versus s" << '\n';
  printVector(_rho_PI);
  std::cout << "rho_k versus s" << '\n';
  printVector(_rho_k);
  std::cout << "phi cell-average" << '\n';
  printVector(_phi_cell.back());
}

//1d vector printing function
void TransportEquation::printVector(const std::vector<double>& vec) const
{
  for (auto col = vec.begin(); col != vec.end(); ++col)
    std::cout << *col << "  ";
  std::cout << std::endl;
  std::cout << std::endl;
}

void TransportEquation::printVector(const std::vector<int>& vec) const
{
  for (auto col = vec.begin(); col != vec.end(); ++col)
    std::cout << *col << "  ";
  std::cout << std::endl;
  std::cout << std::endl;
}

//2d vector printing function
void TransportEquation::printVector(const std::vector<std::vector<double>>& vec) const
{
  for (auto row = vec.begin(); row != vec.end(); ++row)
  {
    for (auto col = row->begin(); col != row->end(); ++col)
       std::cout << *col << "  ";
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

//subtraction of 2 vectors
std::vector<double> TransportEquation::vectorSubtract(std::vector<double> v1, const std::vector<double>& v2) const
{
  if(v1.size() != v2.size()) throw std::invalid_argument( "vectors length not equal!" );
  for(auto i=0; i<v1.size(); ++i) v1[i] -= v2[i];
  return v1;
}

std::vector<std::vector<double>> TransportEquation::vectorSubtract(std::vector<std::vector<double>> v1, const std::vector<std::vector<double>>& v2)
const {
  if(v1.size() != v2.size()) throw std::invalid_argument( "vectors length not equal!" );
  for(auto i=0; i<v1.size(); ++i) {
    for (auto j = 0; j < v1[i].size(); j++) {
      v1[i][j] -= v2[i][j];
    }
  }
  return v1;
}

//calculate the infinite norm of a vector
double TransportEquation::normInf(std::vector<double> vec) const
{
  for (auto col = vec.begin(); col != vec.end(); ++col) *col = std::abs(*col);
  return (*std::max_element(vec.begin(),vec.end()));
}

double TransportEquation::normInf(std::vector<std::vector<double>> vec) const
{
  double temp_max=std::abs(vec[0][0]);
  for (auto row = vec.begin(); row != vec.end(); ++row){
    for (auto col = row->begin(); col != row->end(); ++col)
       if (std::abs(*col) > temp_max) temp_max = std::abs(*col);
  }
  return temp_max;
}

//normlize a vector _phi
void TransportEquation::scalarFluxNormalize(std::vector<std::vector<double>>& vec)
{
  double temp_sum = 0.0;
  for(auto i=0; i<vec.size(); ++i){
    for (size_t g = 0; g < vec[i].size(); g++) {
      temp_sum += _cell_area[i%_num_of_cell_per_pin]*vec[i][g];
    }
  }
  for(auto i=0; i<vec.size(); ++i){
    for (size_t g = 0; g < vec[i].size(); g++) {
      vec[i][g] /= temp_sum;
    }
  }
}
