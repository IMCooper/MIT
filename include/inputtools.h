// deal.II includes:
#include <deal.II/base/parameter_handler.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

// std includes:
#include <fstream>
#include <iostream>

// My includes:
#include <all_data.h>

using namespace dealii;

#ifndef INPUTTOOLS_H
#define INPUTTOOLS_H
namespace InputTools
{
  template <int dim>
  void read_in_mesh (std::string mesh_name,
                     Triangulation<dim> &triangulation);
  
  class ParameterReader : public Subscriptor
  {
  public:
    ParameterReader(ParameterHandler &paramhandler);
    ~ParameterReader();
    void read_parameters (const std::string);
    void read_and_copy_parameters (const std::string);
    void copy_parameters ();
    
  private:
    void declare_parameters();
    void copy_to_equation_data();
    void get_matrix_from_list(std::string entry, FullMatrix<double> &matrix_out, unsigned int matrix_size);
    void get_vector_from_list(std::string entry, Vector<double> &vector_out, unsigned int vector_length);
    ParameterHandler &prm;
  };
  
}
#endif
