// deal.II includes:
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

// std includes:

// My includes:
//#include <backgroundfield.h>
#include <curlfunction.h>

using namespace dealii;

#ifndef MYVECTORTOOLS_H
#define MYVECTORTOOLS_H
namespace MyVectorTools
{
  // Calculates the Hcurl error for the input solution against
  // the input curlFunction<dim> exact solution.
  //
  // NOTE: the solution must correspond to a real/imaginary part of an
  //       FE_Nedelec finite element.
  template<int dim, class DH>
  double calcErrorHcurlNorm(const DH &dof_handler,
                            const Vector<double> &solution,
                            const curlFunction<dim> &exact_solution);
  
  // Functions to return the gradient values of an FE solution at a point:
  // Note: only works for FE_Nedelec elements with 2 vector valued blocks
  //       which correspond to the real & imaginary part.
  template <int dim, class DH>
  void point_gradient(const DH &dof_handler,
                      const Vector<double> &solution,
                      const Point<dim> &point,
                      std::vector<Tensor<1,dim>> &gradients);
  
  // Functions to return the curls of an FE solution at a point:
  // Note: only works for FE_Nedelec elements with 2 vector valued blocks
  //       which correspond to the real & imaginary part.
  template <int dim, class DH>
  void point_curl (const DH &dof_handler,
                   const Vector<double> &solution,
                   const Point<dim> &point,
                   Vector<double> &curl);
}
#endif