// deal.II includes:
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>

// std includes:
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// My includes:
#include <all_data.h>
#include <curlfunction.h>
#include <myvectortools.h>

using namespace dealii;

#ifndef OUTPUTTOOLS_H
#define OUTPUTTOOLS_H
namespace OutputTools
{
  // Postprocessor for use by output_to_vtk
  template <int dim>
  class Postprocessor : public DataPostprocessor<dim>
  {
  public:
    // TODO: add functionality to pass a function class to this.
    Postprocessor (const curlFunction<dim> &exact_solution);
    virtual
    void
    compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                       const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                       const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                       const std::vector<Point<dim> >                  &normals,
                                       const std::vector<Point<dim> >                  &evaluation_points,
                                       std::vector<Vector<double> >                    &computed_quantities) const;
                                       
    virtual std::vector<std::string> get_names () const;
    virtual
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation () const;
    virtual UpdateFlags get_needed_update_flags () const;
    
  private:
    const SmartPointer< const curlFunction<dim> > exact_solution;
    
  };
  
  template<int dim, class DH>
  void output_to_vtk (const DH &dof_handler,
                      const Vector<double> &solution,
                      const std::string &vtk_filename,
                      const curlFunction<dim> &exact_solution);
  template<int dim, class DH>
  void output_to_vtk (const Mapping<dim> &mapping,
                      const DH &dof_handler,
                      const Vector<double> &solution,
                      const std::string &vtk_filename,
                      const curlFunction<dim> &exact_solution);
  
  // TODO ADD output_radial_perturbed_values.
  // Then remove the perturbed parts from the output_radial_values routines.
  /*
  template <int dim, class DH>
  void output_radial_perturbed_values(const DH &dof_handler,
                            const Vector<double> &solution,
                            const perturbedFunction<dim> &boundary_conditions,
                            const Vector<double> &uniform_field,
                            const Point<dim> &end_point,
                            const std::string &filename);*/
  
  template <int dim, class DH>
  void output_radial_values(const DH &dof_handler,
                            const Vector<double> &solution,
                            const perturbedFunction<dim> &boundary_conditions,
                            const Vector<double> &uniform_field,
                            const Point<dim> &end_point,
                            const std::string &filename);
  template <int dim, class DH>
  void output_radial_values(const Mapping<dim> &mapping,
                            const DH &dof_handler,
                            const Vector<double> &solution,
                            const perturbedFunction<dim> &boundary_conditions,
                            const Vector<double> &uniform_field,
                            const Point<dim> &end_point,
                            const std::string &filename);
}
#endif