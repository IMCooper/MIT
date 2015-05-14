#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>


using namespace dealii;

#ifndef CURLFUNCTION_H
#define CURLFUNCTION_H
template<int dim>
class curlFunction : public Function<dim>
{
public:
//   curlFunction ();//unsigned int n_components = dim+dim);
  virtual void curl_value_list (const std::vector<Point<dim> > &points,
                                std::vector<Vector<double> >   &values) const = 0;
                                
  //virtual void perturbed_field_value_list (const std::vector<Point<dim> > &points,
  //                                         std::vector<Vector<double> >   &values) const = 0;
};
#endif

#ifndef PERTURBEDFUNCTION_H
#define PERTURBEDFUNCTION_H
template<int dim>
class perturbedFunction : public curlFunction<dim>
{
public:
//   perturbedFunction ();//unsigned int n_components = dim+dim);
  virtual void curl_value_list (const std::vector<Point<dim> > &points,
                                std::vector<Vector<double> >   &values) const = 0;
                                
  virtual void perturbed_field_value_list (const std::vector<Point<dim> > &points,
                                           std::vector<Vector<double> >   &values) const = 0;
};
#endif