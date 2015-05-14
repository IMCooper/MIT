// deal.II includes:
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

// std includes:
#include <complex>


// My includes:
#include <curlfunction.h>
#include <all_data.h>

using namespace dealii;

#ifndef BACKGROUNDFIELD_H
#define BACKGROUNDFIELD_H
namespace backgroundField
{

  // Dipole field centred at a point
  template<int dim>
  class DipoleSource : public curlFunction<dim>
  {
  public:
    DipoleSource(const Point<dim> &input_source_point,
                 const Tensor<1, dim> &input_coil_direction);
        
    void vector_value_list (const std::vector<Point<dim> > &points,
                            std::vector<Vector<double> >   &values) const;
        
    void curl_value_list (const std::vector<Point<dim> > &points,
                          std::vector<Vector<double> >   &values) const;
                                        
  private:
    Point<dim> source_point;
    Tensor<1, dim> coil_direction;
  };
  
  // Field for a conducting sphere in a uniform background field.
  // includes a function to compute the perturbed field.
  template<int dim>
  class conductingSphere : public perturbedFunction<dim>
  {
  public:
    conductingSphere(double sphere_radius,
                     const Vector<double> &uniform_field);
        
    void vector_value_list (const std::vector<Point<dim> > &points,
                            std::vector<Vector<double> >   &value_list) const;
        
    void curl_value_list (const std::vector<Point<dim> > &points,
                          std::vector<Vector<double> >   &value_list) const;
                          
    void perturbed_field_value_list (const std::vector< Point<dim> > &points,
                                     std::vector<Vector<double> > &value_list) const;
    
    void check_spherical_coordinates(const std::vector< Point<dim> > &points,
                                     std::vector<Vector<double> > &value_list) const;
      
                                        
  private:
    double sphere_radius;
    double constant_p;
    std::complex<double> constant_C;
    std::complex<double> constant_D;
    double constant_B_magnitude;
    Vector<double> uniform_field_re;
    Vector<double> uniform_field_im;    
  };
  // Field for a conducting cube in a uniform background field.
  // includes a function to compute the perturbed field.
  // This is done using the pertubation tensor which must be read in.
  template<int dim>
  class conductingCube : public perturbedFunction<dim>
  {
  public:
    conductingCube(const Vector<double> &uniform_field);
        
    void vector_value_list (const std::vector<Point<dim> > &points,
                            std::vector<Vector<double> >   &value_list) const;
        
    void curl_value_list (const std::vector<Point<dim> > &points,
                          std::vector<Vector<double> >   &value_list) const;
                          
    void perturbed_field_value_list (const std::vector< Point<dim> > &points,
                                     std::vector<Vector<double> > &value_list) const;
      
                                        
  private:
    double cube_size=0.005;
    Vector<double> uniform_field_re;
    Vector<double> uniform_field_im;    
  };
  // Wave propagation
  // E = p*exp(i*k*x), x in R^3
  //
  // with p orthogonal to k, omega = |k| & |p|=1.
  template<int dim>
  class WavePropagation : public curlFunction<dim> 
  {
  public:
    WavePropagation(Vector<double> &k_wave,
                    Vector<double> &p_wave);
    
    void vector_value_list (const std::vector<Point<dim> > &points,
                            std::vector<Vector<double> >   &values) const;
                            

    void curl_value_list(const std::vector<Point<dim> > &points,
                         std::vector<Vector<double> > &value_list) const;
    
  private:  
    Vector<double> k_wave;
    Vector<double> p_wave;
  };
  
  // Function for testing only.
  // Solution is:
  // A = (x^2, y^2, z^2).
  // curlA = (0, 0, 0).
  template<int dim>
  class polynomialTest : public curlFunction<dim>
  {
  public:
    
    polynomialTest();
    
    void vector_value_list (const std::vector<Point<dim> > &points,
                            std::vector<Vector<double> >   &values) const;
                            

    void curl_value_list (const std::vector<Point<dim> > &points,
                          std::vector<Vector<double> > &value_list) const;
  private:  
    Vector<double> k_wave;
    Vector<double> p_wave;
  };
}
#endif