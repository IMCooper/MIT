/*
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>
*/
#include <curlfunction.h>

using namespace dealii;

/*
template<int dim>
class curlFunction : public Function<dim>
{
public:
  virtual void curl_value_list (const std::vector<Point<dim> > &points,
                                std::vector<Vector<double> >   &values) const = 0;
  
  //virtual void perturbed_field_value_list (const std::vector<Point<dim> > &points,
  //                                         std::vector<Vector<double> >   &values) const = 0;
};
*/
// template<int dim>
// curlFunction<dim>::curlFunction()//unsigned int n_components)
// :
// // Function<dim> (n_components)
// {
// }
template class curlFunction<3>;
/*
template<int dim>
class perturbedFunction : public curlFunction<dim>
{
public:
  virtual void curl_value_list (const std::vector<Point<dim> > &points,
                                std::vector<Vector<double> >   &values) const = 0;
                                
  virtual void perturbed_field_value_list (const std::vector<Point<dim> > &points,
                                           std::vector<Vector<double> >   &values) const = 0;
};
*/
// template<int dim>
// perturbedFunction<dim> ::perturbedFunction()//unsigned int n_components)
// :
// // curlFunction<dim> (n_components)
// {}
template class perturbedFunction<3>;