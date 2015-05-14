// deal.II includes:
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

// std includes:

// My includes:

using namespace dealii;
  
#ifndef INVERSESOLVER_H
#define INVERSESOLVER_H
namespace InverseSolver
{
  // Namespace for the inverse solver.
  // Will need to pass a vector of solution vectors
  // and a dof_handler associated with the FE used to
  // compute these solutions.
  // 
  // Note that the solutions will correspond to different
  // BCs but the underlying elements, grid and material
  // parameter distribution will be the same for all of them.
  
  template <int dim>
  void process_mesh(Triangulation<dim> &triangulation);
  
  void initialise_material_parameters();
  
  void setup_coils();
  
  template <int dim>
  void initialise_inverse_problem(Triangulation<dim> &triangulation);
  
  template <int dim, class DH>
  class InverseSolver
  {
  public:
    InverseSolver(std::vector<Vector<double>> &solutions_in,
                  FiniteElement<dim> &fe);
    
    void assemble_sensitivity_matrix(const DH &dof_handler);
    void assemble_sensitivity_rhs(const DH &dof_handler);
    void assemble_regularisation(const DH &dof_handler);
    
    void gauss_newton_solve(const Vector<double> &last_solution,
                            Vector<double> &solution_out);
    
  private:
    const SmartPointer< const FiniteElement<dim> > fe;
    
    std::vector<Vector<double>> solutions;

    
    double regularisation_parameter = InverseProblemData::gn_regularisation_parameter;
    
    
    FullMatrix<double> regularisation_matrix_re;
    FullMatrix<double> regularisation_matrix_im;
    
    FullMatrix<double> sensitivity_matrix_re;
    FullMatrix<double> sensitivity_matrix_im;
    
    Vector<double> sensitivity_rhs_re;
    Vector<double> sensitivity_rhs_im;
  };
}
#endif