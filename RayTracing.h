#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "SegmentsData.h"
#include "PinCellMaterial.h"
#include "LegendrePolynomial.h"

class RayTracing
{
public:
  RayTracing(const int& num_of_lattice_x, const int& num_of_lattice_y) :
  _num_of_lattice_x(num_of_lattice_x), _num_of_lattice_y(num_of_lattice_y)
  {
    _radius = {0.135,0.27,0.405,0.54,0.585};
    _lattice_width = 1.26;
    _lattice_height = 1.26;
    _num_of_lattice_total = _num_of_lattice_x*_num_of_lattice_y;
    _geometry_width = _num_of_lattice_x * _lattice_width;
    _geometry_height = _num_of_lattice_y * _lattice_height;
    _num_of_cells = _radius.size()+1;
    _epsilon_perturb = 1E-12;
    _t_s_original = 0.025;
    setupDirections(16,6);
  }
  SegmentsData trackSegmentation() const;
private:
  std::vector<double> _radius;
  double _lattice_width, _lattice_height, _geometry_width, _geometry_height, _epsilon_perturb, _t_s_original;
  int _num_of_cells, _num_of_lattice_x, _num_of_lattice_y, _num_of_lattice_total;
  std::vector<double> _azim_phi, _polar_theta, _weight_m, _weight_p, _t_s_eff;
  std::vector<int> _n_x, _n_y, _num_of_lines;
  void setupDirections(const unsigned int& M, const unsigned int& P);
  std::vector<double> findUniverse(const std::vector<double>&, std::vector<int>&) const;
  int findCell(const std::vector<double>&) const;
  double computeMinDistance(const std::vector<double>& r_global, const std::vector<double>& direction) const;
  double singleCircleDistance(const std::vector<double>& r_local, const std::vector<double>& direction, const double& r) const;
  double outerRectangleDistance(const std::vector<double>& r_local, const std::vector<double>& direction) const;
  std::vector<double> findIntersectionPoint(const std::vector<double>& r_global, const std::vector<double>& direction) const;
  std::vector<double> getStartingPoint(const int& m, const int& n) const;
  void perturbPoint(std::vector<double>& r_global, const std::vector<double>& direction) const;
};

void RayTracing::setupDirections(const unsigned int& M, const unsigned int& P){
  _azim_phi.resize(M,0); _n_x.resize(M,0); _n_y.resize(M,0); _num_of_lines.resize(M,0); _t_s_eff.resize(M,0);
  for (size_t m = 1; m <= M; m++) { // set up _azim_phi
    _azim_phi[m-1] = 2*M_PI/M*(m-0.5);
    _n_x[m-1] = floor(_geometry_width/_t_s_original*std::abs(sin(_azim_phi[m-1])) )+ 1;
    _n_y[m-1] = floor(_geometry_height/_t_s_original*std::abs(cos(_azim_phi[m-1])) )+ 1;
    _num_of_lines[m-1] = _n_x[m-1] + _n_y[m-1];
    if(_azim_phi[m-1]>0 && _azim_phi[m-1]<M_PI/2)
      _azim_phi[m-1] = atan( _geometry_height*_n_x[m-1]/(_geometry_width*_n_y[m-1]) );
    else if(_azim_phi[m-1]>M_PI/2 && _azim_phi[m-1]<M_PI)
      _azim_phi[m-1] = M_PI - atan( _geometry_height*_n_x[m-1]/(_geometry_width*_n_y[m-1]) );
    else if(_azim_phi[m-1]>M_PI && _azim_phi[m-1]<M_PI*3./2.)
      _azim_phi[m-1] = M_PI + atan( _geometry_height*_n_x[m-1]/(_geometry_width*_n_y[m-1]) );
    else if(_azim_phi[m-1]>M_PI*3./2. && _azim_phi[m-1]<M_PI*2)
      _azim_phi[m-1] = 2*M_PI - atan( _geometry_height*_n_x[m-1]/(_geometry_width*_n_y[m-1]) );
    _t_s_eff[m-1] = _geometry_width/_n_x[m-1] * std::abs(sin(_azim_phi[m-1]));
  }
  _weight_m.resize(M,0);
  for (size_t m = 0; m < M; m++) { // set up _weight_m
    if (m==0) _weight_m[m] = 1./(2*M_PI)*( (_azim_phi[m+1]-_azim_phi[m])/2. + _azim_phi[m] );
    else if(m==M-1) _weight_m[m] = 1./(2*M_PI)*( 2*M_PI - _azim_phi[m] + (_azim_phi[m]-_azim_phi[m-1])/2. );
    else _weight_m[m] = 1./(2*M_PI)*( (_azim_phi[m+1]-_azim_phi[m])/2. + (_azim_phi[m]-_azim_phi[m-1])/2. );

    _weight_m[m] *= 2*M_PI; // to make sure that the weights sum up to 2*pi
  }
  _polar_theta.resize(P,0); _weight_p.resize(P,0);
  LegendrePolynomial LP_obj(-1, 1., P);
  std::vector<double> polar_cos_theta = LP_obj.getRoot();
  for (size_t m = 0; m < polar_cos_theta.size(); m++) {
    _polar_theta[m] = acos(polar_cos_theta[m]);
  }
  _weight_p = LP_obj.getWeight();
}

std::vector<double> RayTracing::findUniverse(const std::vector<double>& r_global, std::vector<int>& current_lattice) const{
  current_lattice.clear();
  current_lattice.resize(2,0);
  current_lattice[0] = floor(r_global[0]/_lattice_width);
  current_lattice[1] = floor(r_global[1]/_lattice_height);
  std::vector<double> r_local = {r_global[0]-(current_lattice[0]+0.5)*_lattice_width, r_global[1]-(current_lattice[1]+0.5)*_lattice_height};
  return r_local;
}

int RayTracing::findCell(const std::vector<double>& r_global) const{
  std::vector<int> current_lattice;
  std::vector<double> r_local = findUniverse(r_global, current_lattice);
  double distance_to_center = sqrt(std::pow(r_local[0],2) + std::pow(r_local[1],2));
  int current_cell;
  if(distance_to_center > *std::max_element(_radius.begin(),_radius.end())){
    current_cell = _radius.size();
  }
  else{
    for (size_t i = 0; i < _radius.size(); i++) {
      if(distance_to_center<_radius[i]){
        current_cell=i;
        break;
      }
    }
  }
  return current_cell;
}

double RayTracing::computeMinDistance(const std::vector<double>& r_global, const std::vector<double>& direction) const{
  double d_min = _lattice_width + _lattice_height;
  std::vector<int> current_lattice;
  std::vector<double> r_local = findUniverse(r_global, current_lattice);
  int current_cell = findCell(r_global);
  std::vector<double> possible_distances;
  if(current_cell==0) return(singleCircleDistance(r_local, direction, _radius[0]));
  else if(current_cell==_radius.size()){
    double d1 = singleCircleDistance(r_local, direction, _radius[current_cell-1]);
    double d2 = outerRectangleDistance(r_local, direction);
    return (d1<d2 ? d1:d2);
  }
  else{
    double d1 = singleCircleDistance(r_local, direction, _radius[current_cell-1]);
    double d2 = singleCircleDistance(r_local, direction, _radius[current_cell]);
    return (d1<d2 ? d1:d2);
  }
}

double RayTracing::outerRectangleDistance(const std::vector<double>& r_local, const std::vector<double>& direction) const{
  std::vector<double> line_segments;
  double temp_line;
  temp_line = sqrt(1+std::pow(direction[1]/direction[0], 2))*(_lattice_width/2-r_local[0]);
  if(temp_line>10*_epsilon_perturb) line_segments.push_back(temp_line);
  temp_line = sqrt(1+std::pow(direction[1]/direction[0], 2))*(_lattice_width/2+r_local[0]);
  if(temp_line>10*_epsilon_perturb) line_segments.push_back(temp_line);
  temp_line = sqrt(1+std::pow(direction[0]/direction[1], 2))*(_lattice_height/2-r_local[1]);
  if(temp_line>10*_epsilon_perturb) line_segments.push_back(temp_line);
  temp_line = sqrt(1+std::pow(direction[0]/direction[1], 2))*(_lattice_height/2+r_local[1]);
  if(temp_line>10*_epsilon_perturb) line_segments.push_back(temp_line);
  double distance;
  for (size_t i = 0; i < line_segments.size(); i++) {
    std::vector<double> next_point = {r_local[0]+line_segments[i]*direction[0], r_local[1]+line_segments[i]*direction[1]};
    if( (std::abs(next_point[0]) <= _lattice_width/2+10*_epsilon_perturb) && (std::abs(next_point[1]) <= _lattice_height/2+10*_epsilon_perturb))
      if( (std::abs(std::abs(next_point[0])-_lattice_width/2) < _epsilon_perturb) ||  (std::abs(std::abs(next_point[1])-_lattice_height/2) < 1E-2*_epsilon_perturb))
        return line_segments[i];
  }
  return (_lattice_width+_lattice_height);
}

double RayTracing::singleCircleDistance(const std::vector<double>& r_local, const std::vector<double>& direction, const double& r) const{
  auto d_sqr_lambda = [](const auto& a) { return (std::pow(a[0],2) + std::pow(a[1],2));};
  double d_line_to_center = std::abs(r_local[1]*direction[0]-r_local[0]*direction[1])/d_sqr_lambda(direction);
  double d_point_to_center = sqrt(d_sqr_lambda(r_local));
  if(d_line_to_center < r){
    double half_chord = sqrt( std::pow(r,2) - std::pow(d_line_to_center,2) );
    double temp = sqrt( std::pow(d_point_to_center,2) - std::pow(d_line_to_center,2) );
    std::vector<double> line_segments;
    if(d_point_to_center < r){
      line_segments.push_back(half_chord+temp);
      if(std::abs(half_chord-temp) > 10*_epsilon_perturb)
        line_segments.push_back(std::abs(half_chord-temp));
    }
    else{
      line_segments.push_back(half_chord+temp);
      if(std::abs(half_chord-temp) > 10*_epsilon_perturb)
        line_segments.push_back(std::abs(temp-half_chord));
    }
    //check if on the right direction
    double min_distance=_lattice_width+_lattice_height;
    for (size_t i = 0; i < line_segments.size(); i++) {
      std::vector<double> next_point = {r_local[0]+line_segments[i]*direction[0], r_local[1]+line_segments[i]*direction[1]};
      if(std::abs(sqrt(d_sqr_lambda(next_point))-r) <_epsilon_perturb){
        if(line_segments[i]<min_distance) min_distance=line_segments[i];
      }
    }
    return min_distance;
  }
  else{
    return (_lattice_width+_lattice_height);//just a very big number;
  }
}

std::vector<double> RayTracing::findIntersectionPoint(const std::vector<double>& r_global, const std::vector<double>& direction) const{
  double d = computeMinDistance(r_global, direction);
  double next_x=r_global[0]+d*direction[0], next_y=r_global[1]+d*direction[1];
  if(std::abs(next_x)<_epsilon_perturb*1E-3) next_x=0;
  if(std::abs(next_y)<_epsilon_perturb*1E-3) next_y=0;
  std::vector<double> temp_vec({next_x, next_y});
  return temp_vec;
}

std::vector<double> RayTracing::getStartingPoint(const int& m, const int& n) const{
  std::vector<double> startPoint;
  if(_azim_phi[m]>0 && _azim_phi[m]<M_PI/2){
    if(n<_n_x[m]) startPoint = { _t_s_eff[m]/std::abs(sin(_azim_phi[m])) * (0.5+n), 0 };
    else if(n<_n_x[m]+_n_y[m]) startPoint = { 0, _t_s_eff[m]/std::abs(cos(_azim_phi[m])) * (0.5+n-_n_x[m]) };
    else throw std::invalid_argument( "Invalid number of n!" );
  }
  else if(_azim_phi[m]>M_PI/2 && _azim_phi[m]<M_PI){
    if(n<_n_x[m]) startPoint = { _geometry_width - _t_s_eff[m]/std::abs(sin(_azim_phi[m])) * (0.5+n), 0 };
    else if(n<_n_x[m]+_n_y[m]) startPoint = { _geometry_width, _t_s_eff[m]/std::abs(cos(_azim_phi[m])) * (0.5+n-_n_x[m]) };
    else throw std::invalid_argument( "Invalid number of n!" );
  }
  else if(_azim_phi[m]>M_PI && _azim_phi[m]<M_PI*3./2.){
    if(n<_n_x[m]) startPoint = { _geometry_width - _t_s_eff[m]/std::abs(sin(_azim_phi[m])) * (0.5+n), _geometry_height };
    else if(n<_n_x[m]+_n_y[m]) startPoint = { _geometry_width, _geometry_height - _t_s_eff[m]/std::abs(cos(_azim_phi[m])) * (0.5+n-_n_x[m]) };
    else throw std::invalid_argument( "Invalid number of n!" );
  }
  else if(_azim_phi[m]>M_PI*3./2. && _azim_phi[m]<M_PI*2){
    if(n<_n_x[m]) startPoint = { _t_s_eff[m]/std::abs(sin(_azim_phi[m])) * (0.5+n), _geometry_height };
    else if(n<_n_x[m]+_n_y[m]) startPoint = { 0, _geometry_height - _t_s_eff[m]/std::abs(cos(_azim_phi[m])) * (0.5+n-_n_x[m]) };
    else throw std::invalid_argument( "Invalid number of n!" );
  }
  else throw std::invalid_argument( "Invalid value of _azim_phi[m]!" );
  return startPoint;
}

SegmentsData RayTracing::trackSegmentation() const{
  SegmentsData SD_obj(_azim_phi.size(), _num_of_lines);
  PinCellMaterial PCM_obj(_num_of_lattice_x, _num_of_lattice_y);
  SD_obj._azim_phi = _azim_phi;
  SD_obj._polar_theta = _polar_theta;
  SD_obj._weight_az = _weight_m;
  SD_obj._weight_po = _weight_p;
  SD_obj._cell_to_material.resize(_num_of_cells*_num_of_lattice_total);
  for (auto n_y = 0; n_y < _num_of_lattice_y; n_y++) {
    for (auto n_x = 0; n_x < _num_of_lattice_x; n_x++) {
      for (size_t p = 0; p < _num_of_cells; p++) {
        std::vector<int> current_lattice = {n_x,n_y};
        std::vector<double> point_in_cell;
        if(p < _num_of_cells-1) point_in_cell = {_radius[p]-_epsilon_perturb, 0};
        else point_in_cell = {_radius[p-1]+_epsilon_perturb, 0};
        SD_obj._cell_to_material[ p + _num_of_cells*(n_x+n_y*_num_of_lattice_x) ] = PCM_obj.findFlatSource(point_in_cell, current_lattice);
      }
    }
  }

  for (size_t m = 0; m < _azim_phi.size(); m++) {
    for (size_t n = 0; n < _n_x[m]+_n_y[m]; n++) {
      SD_obj._width[m][n] = _t_s_eff[m];
      std::vector<double> current_point = getStartingPoint(m,n);
      SD_obj._point[m][n].push_back(current_point);
      std::vector<double> direction = {cos(_azim_phi[m]), sin(_azim_phi[m])};
      perturbPoint(current_point, direction);
      while (current_point[0]>0 && current_point[0]<_geometry_width && current_point[1]>0 && current_point[1]<_geometry_height) {
        std::vector<double> next_point = findIntersectionPoint(current_point, direction);
        SD_obj._point[m][n].push_back(next_point);
        perturbPoint(next_point, direction);
        std::vector<int> current_lattice;
        std::vector<double> r_local = findUniverse(current_point, current_lattice);
        SD_obj._lattice[m][n].push_back(current_lattice[0] + current_lattice[1]*_num_of_lattice_x);
        SD_obj._material[m][n].push_back(PCM_obj.findFlatSource(r_local, current_lattice));
        SD_obj._cell[m][n].push_back(findCell(current_point)+_num_of_cells*(current_lattice[0] + current_lattice[1]*_num_of_lattice_x));
        SD_obj._length[m][n].push_back(sqrt( std::pow(next_point[0]-current_point[0],2)+std::pow(next_point[1]-current_point[1],2) ));
        current_point = next_point;
      }
    }
  }
  // find corresponding reflected point indexes
  for (size_t m = 0; m < _azim_phi.size(); m++) {
    for (size_t n = 0; n < _n_x[m]+_n_y[m]; n++) {
      std::vector<double> last_point = SD_obj._point[m][n].back();
      std::vector<unsigned long int> possible_index_m;
      if(_azim_phi[m]>0 && _azim_phi[m]<M_PI) possible_index_m = {_azim_phi.size()-m-1, _azim_phi.size()/2-m-1};
      else possible_index_m = {_azim_phi.size()-m-1, _azim_phi.size()*3/2-m-1};
      for (size_t mm = 0; mm < possible_index_m.size(); mm++) {
        for (size_t nn = 0; nn < _n_x[possible_index_m[mm]]+_n_y[possible_index_m[mm]]; nn++) {
          std::vector<double> temp_point = SD_obj._point[possible_index_m[mm]][nn][0];
          if(sqrt(std::pow(temp_point[0]-last_point[0],2) + std::pow(temp_point[1]-last_point[1],2)) < 10*_epsilon_perturb){
            SD_obj._reflect_point_m_az[m][n] = possible_index_m[mm];
            SD_obj._reflect_point_n[m][n] = nn;
            break;
          }
        }
      }
    }
  }
  return SD_obj;
}

void RayTracing::perturbPoint(std::vector<double>& r_global, const std::vector<double>& direction) const{
  r_global = {r_global[0]+direction[0]*_epsilon_perturb, r_global[1]+direction[1]*_epsilon_perturb};;
}

#endif
