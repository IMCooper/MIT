#include <deal.II/base/quadrature.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <curlfunction.h>

using namespace dealii;

#ifndef VOLTAGE_H
#define VOLTAGE_H
namespace Voltage
{
  template <int dim>
  class Voltage
  {
    // Class to calculate the voltage for a particular sensor coil.
    // The sensor coil must be defined by a point (centre of the coil)
    // and a direction (moment), which is a vector [(1,dim) tensor]
    // 
    // The approximation implemented is to assume the coil has zero width
    // and that we may perform a surface integral over the coil:
    //
    // voltage = \int_{coil} E*tangent dS
    //
    // However, our solution is A, with E = -i*omega*A so,
    //
    // voltage = -i*omega*(\int_{coil} A*tangent dS.
  public:
    Voltage(const Point<dim> &input_measurement_point,
            const Tensor<1, dim> &input_coil_direction,
            const unsigned int quad_order_in = 20);
    
    void calculateVoltage(const Vector<double> &solution_in,
                          const DoFHandler<dim> &dof_handler,
                          double &voltage_out_re,
                          double &voltage_out_im);
                                        
  private:
    unsigned int quad_order;
    
    Tensor<1,dim> rotate_coil(const Tensor<1, dim> &input_point,
                              const Tensor<1, dim> &axis,
                              const double &angle);
                           
    
    // Measurement coil data:
    Point<dim> measurement_point;
    Tensor<1, dim> coil_direction;
    
    unsigned int n_quadrature_points;
    std::vector< Point<dim> > coil_quad_points;
    std::vector< Tensor<1,dim> > coil_tangent_vectors;
    
    // Reference coil data:
    std::vector<double> reference_quad_weights;
    Vector<double> reference_quad_angles;
    std::vector< Point<dim> > reference_quad_points;
    std::vector< Tensor<1,dim> > reference_tangent_vectors;
    double reference_jacobian = 2.0*numbers::PI*SensorCoilData::coil_radius;
    Tensor<1,3> reference_coil_direction;
  };  
  
  
  template<int dim, class DH>
  void alternateMeasureVoltage(/*const std::vector<Point<dim>> &all_exciter_centres,
                               const std::vector<Tensor<1, dim>> &all_exciter_directions,*/
                               const std::vector<Point<dim>> &all_sensor_centres,
                               const std::vector<Tensor<1, dim>> &all_sensor_directions,
                               const std::vector<Vector<double>> &all_solutions,
                               const std::vector<SmartPointer<backgroundField::curlFunction<dim>>> &all_boundary_functions,
                               const DH &dof_handler,
                               std::vector<Vector<double>> &measured_voltages_re,
                               std::vector<Vector<double>> &measured_voltages_im);
}
#endif