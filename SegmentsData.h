#ifndef SEGMENTSDATA_H
#define SEGMENTSDATA_H

#include <vector>

class SegmentsData
{
public:
  SegmentsData(){}
  SegmentsData(const int& M, const std::vector<int>& num_of_lines){
    _length.resize(M); _material.resize(M); _cell.resize(M); _point.resize(M); _width.resize(M);
    _reflect_point_m_az.resize(M); _reflect_point_n.resize(M); _lattice.resize(M);
    for (size_t m = 0; m < M; m++) {
      _length[m].resize(num_of_lines[m]);
      _material[m].resize(num_of_lines[m]);
      _cell[m].resize(num_of_lines[m]);
      _point[m].resize(num_of_lines[m]);
      _width[m].resize(num_of_lines[m]);
      _reflect_point_m_az[m].resize(num_of_lines[m]);
      _reflect_point_n[m].resize(num_of_lines[m]);
      _lattice[m].resize(num_of_lines[m]);
    }
  }
  std::vector<std::vector<std::vector<double>>> _length;
  std::vector<std::vector<std::vector<int>>> _material, _cell;
  std::vector<std::vector<std::vector<std::vector<double>>>> _point;
  std::vector<std::vector<double>> _width;
  std::vector<std::vector<int>> _reflect_point_m_az, _reflect_point_n;
  std::vector<double> _azim_phi, _polar_theta, _weight_az, _weight_po;
  std::vector<int> _cell_to_material;
  std::vector<std::vector<std::vector<int>>> _lattice;
};

#endif
