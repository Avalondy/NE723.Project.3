#ifndef PINCELLMATERIAL_H
#define PINCELLMATERIAL_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

class PinCellMaterial
{
public:
  PinCellMaterial(){}
  PinCellMaterial(const int& num_of_lattice_x, const int& num_of_lattice_y){
    _material_data_file = "material.Input.txt";
    _num_of_table = 7;
    _G = 7;
    _fuel_radius = 0.54;
    _pin_cell_width = 1.26;
    // readMaterialCrosssections();
    _num_of_pin_cell = num_of_lattice_x*num_of_lattice_y;
    if(_num_of_pin_cell == 1) _pin_cell_type = {{MaterialEnum::UO2}};
    else if(_num_of_pin_cell == 9){
      _pin_cell_type = {{MaterialEnum::MOX43, MaterialEnum::UO2, MaterialEnum::MOX43},
                        {MaterialEnum::UO2, MaterialEnum::GuideTube, MaterialEnum::UO2},
                        {MaterialEnum::MOX87, MaterialEnum::UO2, MaterialEnum::MOX70}};
    }
    else std::cout << "Invalid value of _num_of_pin_cell" << '\n';

    fillUpCrossections();
  }

  std::vector<std::vector<double>> _Sigma_t, _Sigma_f, _nu_f, _chi;
  std::vector< std::vector< std::vector<double> > > _Sigma_s;

  int findFlatSource(const std::vector<double>& r_local, const std::vector<int>& current_lattice) const;
  void fillUpCrossections();
private:
  enum class MaterialEnum {UO2, MOX43, MOX70, MOX87, GuideTube, Moderator};
  int _num_of_material_enum = 6;
  double _fuel_radius, _pin_cell_width;
  std::vector<std::vector<MaterialEnum>> _pin_cell_type;
  std::string _material_data_file;
  int _num_of_table, _num_of_pin_cell, _G;
  std::vector< std::vector< std::vector<double> > > _cross_sections_table, _scattering_table;
  void readMaterialCrosssections();
};

//read the C5G7 tables
void PinCellMaterial::readMaterialCrosssections(){
  _cross_sections_table.resize(_num_of_table);
  _scattering_table.resize(_num_of_table);
  std::ifstream data_file(_material_data_file);
  std::string line;
  if (data_file.is_open()){
    int i_table=-1, i_row=0, i_col=0, which_table;
    while ( getline (data_file,line) ){
      if(line.empty()||line.rfind("Groups#", 0) == 0) continue;
      if(line.rfind("Table", 0) == 0){
        ++i_table;
        which_table=1;
        continue;
      }
      else if(line.rfind("Scattering Block", 0) == 0){
        which_table=2;
        continue;
      }
      std::string dump;
      double number;
      std::vector<double> line_number_array;
      std::istringstream iss(line);
      iss >> dump;
      while(iss >> number) line_number_array.push_back(number);
      if(which_table==1)
        _cross_sections_table[i_table].push_back(line_number_array);
      else if(which_table==2)
        _scattering_table[i_table].push_back(line_number_array);
    }

    data_file.close();
  }
  else std::cout << "Unable to open file" << '\n';
}

int PinCellMaterial::findFlatSource(const std::vector<double>& r_local, const std::vector<int>& current_lattice)const{
  double distance_to_center = sqrt(std::pow(r_local[0],2) + std::pow(r_local[1],2));
  if(distance_to_center <= _fuel_radius){
    return static_cast<int>(_pin_cell_type[ current_lattice[1] ][ current_lattice[0] ]);
  }
  else if(std::abs(r_local[0])<=_pin_cell_width/2 && std::abs(r_local[1])<=_pin_cell_width/2){
    return static_cast<int>(MaterialEnum::Moderator);
  }
  else{
    std::cout << "Invalid value of r_local!" << '\n';
    return -1;
  }
}

void PinCellMaterial::fillUpCrossections(){
  readMaterialCrosssections();
  _Sigma_t.resize(_num_of_material_enum, std::vector<double>(_G, 0));
  _Sigma_f.resize(_num_of_material_enum, std::vector<double>(_G, 0));
  _nu_f.resize(_num_of_material_enum, std::vector<double>(_G, 0));
  _chi.resize(_num_of_material_enum, std::vector<double>(_G, 0));
  _Sigma_s.resize(_num_of_material_enum);
  for (size_t material_type_int = 0; material_type_int < _num_of_material_enum; material_type_int++) {
    MaterialEnum material_type = static_cast<MaterialEnum>(material_type_int);
    int which_table=-1;
    if(material_type==MaterialEnum::UO2) which_table=0;
    else if(material_type==MaterialEnum::MOX43) which_table=1;
    else if(material_type==MaterialEnum::MOX70) which_table=2;
    else if(material_type==MaterialEnum::MOX87) which_table=3;
    else if(material_type==MaterialEnum::GuideTube) which_table=5;
    else if(material_type==MaterialEnum::Moderator) which_table=6;
    else std::cout << "Invalid value of material_type!" << '\n';

    for (auto g = 0; g < _G; g++) {
      _Sigma_t[material_type_int][g] = _cross_sections_table[which_table][g][1];
      if(which_table==5 || which_table==6){
        _Sigma_f[material_type_int][g] = 0.;
           _nu_f[material_type_int][g] = 0.;
            _chi[material_type_int][g] = 0.;
      }
      else{
        _Sigma_f[material_type_int][g] = _cross_sections_table[which_table][g][4];
           _nu_f[material_type_int][g] = _cross_sections_table[which_table][g][5];
            _chi[material_type_int][g] = _cross_sections_table[which_table][g][6];
      }
    }
    _Sigma_s[material_type_int] = _scattering_table[which_table];
  }
}

#endif
