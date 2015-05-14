#include <voltage.h>

  using namespace dealii;

namespace Voltage
{ // Useful functions:
  // 
  // GridTools::find_active_cell_around_point - this finds a cell containing a given point.
  // 
  // VectorTools::point_value - evaluate the solution at a given point, think it fails if you're not in the correct cell already.
  // 
  // useful threads:
  // https://groups.google.com/forum/#!searchin/dealii/evaluate$20solution$20point/dealii/ShYw9VSOz1A/DQUnJ6hHl9IJ
  // https://groups.google.com/forum/#!searchin/dealii/evaluate$20solution$20point/dealii/ZXLuabqcpT4/9W9gfj0yziMJ
  // 

    
  template <int dim>
  Voltage<dim>::Voltage(const Point<dim> &input_measurement_point,
                        const Tensor<1, dim> &input_coil_direction,
                        const unsigned int quad_order_in)
  :
  quad_order(quad_order_in),
  measurement_point(input_measurement_point),
  coil_direction(input_coil_direction)  
  {
    // Setup standard coil centred at (0,0,0) with radius coil_radius, direction straight up (0,0,1)
    // and with points mapped from a QGauss quad rule (0,1) to the disc (-pi,pi)
    
    // Mid-point quad rule (REMOVED):
//     n_quadrature_points = 10;
    
    // Gauss quad rule:
    QGauss<1> reference_quadrature(quad_order);
    n_quadrature_points = reference_quadrature.size();  

    // Map quad points from 0, 1 to -pi, pi.
    std::vector< Point<1> > temp_points (n_quadrature_points);
    temp_points = reference_quadrature.get_points();
    
    // Gauss weights:
    reference_quad_weights.resize(n_quadrature_points);
    reference_quad_weights = reference_quadrature.get_weights();
    
    // Calculate quadrature angle and points
    // and tangent vectors to the disc at these points.:
    
    reference_quad_angles.reinit(n_quadrature_points);
    reference_quad_points.resize(n_quadrature_points);
    reference_tangent_vectors.resize(n_quadrature_points);
    
    for (unsigned int i = 0; i < n_quadrature_points; ++i)
    {
      // Mid point angles:
      // reference_quad_angles(i) = 2.0*numbers::PI/n_quadrature_points*((double)i);
      
      // Mid point weights:
      // reference_quad_weights[i] = 1.0/(double)n_quadrature_points;
      
      // Gauss point angles:
      reference_quad_angles(i) = numbers::PI*(2.0*temp_points[i](0) - 1.0);      
      
      // Points:
      reference_quad_points[i](0) = SensorCoilData::coil_radius*cos(reference_quad_angles(i));
      reference_quad_points[i](1) = SensorCoilData::coil_radius*sin(reference_quad_angles(i));
      reference_quad_points[i](2) = 0.0;
      
      // Tangents::
      reference_tangent_vectors[i][0] = -sin(reference_quad_angles(i));
      reference_tangent_vectors[i][1] =  cos(reference_quad_angles(i));
      reference_tangent_vectors[i][2] =  0.0;      

    }    

    
    // Store reference coil direction:
    reference_coil_direction[0] = 0.0;
    reference_coil_direction[1] = 0.0;
    reference_coil_direction[2] = 1.0;
        
    // Now calculate rotation from the input coil to the standard coil:
    // This can be used to rotate all the quadrature points on the standard coil
    // to the rotation of the required measurement coil which can then be translated
    //    according to the position of the measurement coil.
    //
    // Can do this by knowing the rotation angle and the rotation axis and then
    //    apply Rodrigues rotation formula,
    //    (rotating vector v about an axis k at angle th)
    //     v_rot = v*cos(th) + cross(axis,v)*sin(th( + k*dot(k,v)*(1-cos(th)).
    
    // rotation angle:
    
    Tensor<1, dim> normalised_coil_direction( coil_direction/coil_direction.norm());
    double rotation_angle = acos(normalised_coil_direction*reference_coil_direction);
    
    // rotation axis:
    Tensor<1, dim> rotation_axis;
    cross_product(rotation_axis,
                  reference_coil_direction,
                  normalised_coil_direction);
    rotation_axis = rotation_axis/rotation_axis.norm ();
    
    /*
     * We now need to rotate and translate the quad points
     * and only rotate the tangent vectors:
     */
    coil_quad_points.resize(n_quadrature_points);
    coil_tangent_vectors.resize(n_quadrature_points);
    for (unsigned int i = 0; i < n_quadrature_points; ++i)
    {
      Point<dim> temp_point(rotate_coil(reference_quad_points[i], rotation_axis, rotation_angle));
      coil_quad_points[i] = temp_point + measurement_point;
      
      Tensor<1, dim> temp_tensor(rotate_coil(reference_tangent_vectors[i], rotation_axis, rotation_angle));
      coil_tangent_vectors[i] = temp_tensor;
    }   
  }
  
  // Voltage members: 
  template <int dim>
  Tensor<1, dim> Voltage<dim>::rotate_coil(const Tensor<1, dim> &input_tensor,
                                           const Tensor<1, dim> &axis,
                                           const double &angle)
  {
    // Rotates a vector (point) about an axis by an angle using
    // Rodrigues rotation formula:
    // (rotating vector v about an axis k at angle th)
    //    v_rot = v*cos(th) + cross(axis,v)*sin(th) + k*dot(k,v)*(1-cos(th)).
    
    Tensor<1,dim> cross_prod;
    cross_product(cross_prod, axis, input_tensor);
    Tensor<1,dim> output_tensor(
      input_tensor*cos(angle)
      + cross_prod*sin(angle)
      + axis*(axis*input_tensor)*(1-cos(angle))
    );
    return output_tensor;
  }
  
  template <int dim>
  void Voltage<dim>::calculateVoltage (const Vector<double> &solution_in,
                                       const DoFHandler<dim> &dof_handler,
                                       double &voltage_out_re,
                                       double &voltage_out_im)
  {
    double integral_re = 0.0;
    double integral_im = 0.0;
    std::vector<Vector<double>> solution_at_quad_point (n_quadrature_points,
                                                        Vector<double> (dof_handler.get_fe().n_components()) );
   
    // Recover solution values at each quad point around the coil
    // then compute the integral.
    //
    // The voltages are:
    // real part:      omega*(int_{coil} A_im*tangent dS
    // imaginary part: -omega*(int_{coil} A_re*tangent dS.

    for (unsigned int i = 0; i < n_quadrature_points; ++i)
    {
      VectorTools::point_value (dof_handler,
                                solution_in,
                                coil_quad_points[i],
                                solution_at_quad_point[i]);
      
      Tensor<1,dim> solution_re;
      Tensor<1,dim> solution_im;
      for (unsigned int d = 0; d < dim; ++d)
      {
        solution_re[d] = solution_at_quad_point[i](d);
        solution_im[d] = solution_at_quad_point[i](3+d);
      }
      
      integral_re += solution_re*coil_tangent_vectors[i]*reference_quad_weights[i];
      integral_im += solution_im*coil_tangent_vectors[i]*reference_quad_weights[i];      
      
    }
    integral_re = reference_jacobian*integral_re*EquationData::param_omega;//( SensorCoilData::coil_radius*SensorCoilData::coil_radius*numbers::PI)*2.0;
    integral_im = reference_jacobian*integral_im*EquationData::param_omega;//( SensorCoilData::coil_radius*SensorCoilData::coil_radius*numbers::PI)*2.0;

    // pass to output:
    voltage_out_re = integral_im;
    voltage_out_im = -integral_re;
  }
  template class Voltage<3>;
  
  template<int dim, class DH>
  void alternateVoltage(/*const std::vector<Point<dim>> &all_exciter_centres,
                        const std::vector<Tensor<1, dim>> &all_exciter_directions,*/
                        const std::vector<Point<dim>> &all_sensor_centres,
                        const std::vector<Tensor<1, dim>> &all_sensor_directions,
                        const std::vector<Vector<double>> &all_solutions,
                        const std::vector<SmartPointer<backgroundField::curlFunction<dim>>> &all_boundary_functions,
                        const DH &dof_handler,
                        std::vector<Vector<double>> &measured_voltages_re,
                        std::vector<Vector<double>> &measured_voltages_im)
  {
    // Routine to caluculate the voltage using an alternative method to
    // the surface contour integral.
    // 
    // To calculate the voltage in coil i centred at the point x
    // for the excitation coil j, we use:
    
    // V|_{i} (x) = m_{i}*( H^{j}_{0}(x) + H_{\Delta}(x) )
    //
    // - H^{j}_{0} = D^{2}G(x,z)*m_{e} is the magnetic field for the exciter coil.
    // - m_{i} is the dipole moment of the measurement coil i.
    // - H_{\Delta}(x) is the magnetic field due to the presence of the conductor
    // - the conducting region is denoted \Omega_{c}
    //
    // H_{\Delta} (x) = -1/(4*pi)*\int_{\Omega_{c}} (x-y)/(|x-y|^{3})*J(y)dy
    //
    // here J(y) = sigma*E=-i*omega*sigma*A. (i.e. the eddy current).
    //
    // i.e. x-y is the radial vector from the sensor coil.
    //      Since the outer region lies in a region with sigma=0, outside \Omega_{c}.
    
    // Check that all the vector input make sense:
    // Let S=number of sensor coils
    //     E=number of excitation coils.
    //
    // we want:
    // all_sensor_centres & all_sensor_directions to be of length S
    // all_solutions & all_boundary_functions to be of length E
    //
    // We resize the output vectors accordingly.
    /*
    Assert(all_exciter_centres.size() == all_exciter_directions.size(),
           ExcDimensionMismatch(all_exciter_centres.size(),
                                all_exciter_directions.size()));
    Assert(all_solutions.size() == all_exciter_centres.size(),
           ExcDimensionMismatch(all_solutions.size(),
                                all_exciter_centres.size()));
    */
    Assert(all_sensor_centres.size() == all_sensor_directions.size(),
           ExcDimensionMismatch(all_sensor_centres.size(),
                                all_sensor_directions.size()));
    Assert(all_solutions.size() == all_boundary_functions.size(),
           ExcDimensionMismatch(all_solutions.size(),
                                all_boundary_functions.size()));
    
    
    
    
    // store the number of excitation/sensor coils:
    unsigned int n_sensors = all_sensor_centres.size();
    unsigned int n_exciters = all_solutions.size();
    
    // resize the output vectors:
    measured_voltages_re.resize(n_exciters);
    measured_voltages_im.resize(n_exciters);
    for (unsigned int exciters=0; exciters<n_exciters; ++exciters)
    {
      measured_voltages_re[exciters].reinit(n_sensors);
      measured_voltages_im[exciters].reinit(n_sensors);
    }
    
    // First generate H_0 for each excitation:
    std::vector<std::vector<Vector<double>>> all_h0_at_centres(n_exciters, std::vector<Vector<double>> (n_sensors,Vector<double> (dim+dim)));
    for (unsigned int exciter=0; exciter<n_exciters; ++exciter)
    {
      // Using vector of backgroundField::curlFunction objects:
      all_boundary_functions[exciter]->curl_value_list(all_sensor_centres,
                                                all_h0_at_centres[exciter]);
      /*
      backgroundField::DipoleSource<dim> boundary_condition(all_exciter_centres[exciter],
                                                            all_exciter_directions[exciter]);
      boundary_condition.curl_value_list(all_sensor_centres,
                                         all_h0_at_centres[exciter]);
      */
      
    }
    // Now calculate each Hdelta:
    // First setup the finite element:
    const FiniteElement<dim> &fe = dof_handler.get_fe ();
    const unsigned int quad_order = 2*fe.degree + 1;    
    QGauss<dim>  quadrature_formula(quad_order);
    const unsigned int n_q_points = quadrature_formula.size();
    FEValues<dim> fe_values(fe, quadrature_formula,
                            update_values | update_quadrature_points | update_JxW_values); // gradients not needed.
    
    const FEValuesExtractors::Vector E_re(0);
    const FEValuesExtractors::Vector E_im(dim);
    
    std::vector<std::vector<Tensor<1,dim> > > all_hdelta_re(n_exciters,
                                              std::vector<Tensor<1,dim>> (n_sensors));
    std::vector<std::vector<Tensor<1,dim> > > all_hdelta_im(n_exciters,
                                              std::vector<Tensor<1,dim>> (n_sensors));
    
    for (unsigned int exciter=0; exciter<n_exciters; ++exciter)
    {
      for (unsigned int sensor=0; sensor<n_sensors; ++sensor)
      {
        all_hdelta_re[exciter][sensor]=0;
        all_hdelta_im[exciter][sensor]=0;
      }
    }
    
    double temp_integral=0.0; // check TESTING
    Point<dim> temp_point(0.0,0.0,0.0);
    
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      // Check if in conducting cell
      if (EquationData::param_sigma(cell->material_id()) > EquationData::param_regularisation)
      {
        double current_sigma = EquationData::param_sigma(cell->material_id());
        fe_values.reinit (cell);
        // Get the value of J for each excitation at each cell quad point
        // Note, J = sigma*E
        // but we have A, E=-i*omega*A
        // so J = -i*omega*sigma*A.
        std::vector<std::vector<Tensor<1,dim> > > cell_solutions_re(n_exciters,
                                                                  std::vector<Tensor<1,dim>> (n_q_points));
        std::vector<std::vector<Tensor<1,dim> > > cell_solutions_im(n_exciters,
                                                                  std::vector<Tensor<1,dim>> (n_q_points));
        for (unsigned int exciter=0; exciter < n_exciters; ++exciter)
        {
          // Calculate real & imaginary parts of E (note switch of A im & re to get E here).
          fe_values[E_re].get_function_values(all_solutions[exciter],
                                              cell_solutions_im[exciter]);
          fe_values[E_im].get_function_values(all_solutions[exciter],
                                              cell_solutions_re[exciter]);
        }
        // Store all quadrature points:        
        std::vector<Point<dim>> all_quad_points = fe_values.get_quadrature_points();
        
        // Loop over all quad points & calculate the contribution to each hdelta
        double local_int=0;
        for (unsigned int q=0; q<n_q_points; ++q)
        {
          // TESTING:
          Point<dim> cyl_temp (all_quad_points[q](0), all_quad_points[q](1), 0.0);
          double temp_radius = cyl_temp.distance(temp_point);
          local_int += temp_radius*fe_values.JxW(q);
          
          
          // calculate (x-y)/(|x-y|^{3}):
          std::vector<Tensor<1,dim>> scaled_vectors(n_sensors);
          for (unsigned int sensor=0; sensor<n_sensors; ++sensor)
          {
            Tensor<1,dim> shifted_point = all_sensor_centres[sensor] - all_quad_points[q];
            double rad = all_sensor_centres[sensor].distance(all_quad_points[q]);
            scaled_vectors[sensor]=shifted_point/(rad*rad*rad);
            for (unsigned int exciter=0; exciter<n_exciters; ++exciter)
            {
              Tensor<1,dim> cross_re;
              Tensor<1,dim> cross_im;
              cross_product(cross_re,
                            scaled_vectors[sensor],
                            cell_solutions_re[exciter][q]);
              cross_product(cross_im,
                            scaled_vectors[sensor],
                            cell_solutions_im[exciter][q]);
              all_hdelta_re[exciter][sensor] += current_sigma*cross_re*fe_values.JxW(q);
              all_hdelta_im[exciter][sensor] += -current_sigma*cross_im*fe_values.JxW(q);
            }
          }
        }
        temp_integral +=local_int;
      }
    }
    std::cout << "Computed integral: " << temp_integral << std::endl;
    // Now calculate voltage
    double integral_factor = -EquationData::param_omega/(4.0*numbers::PI);
    for (unsigned int exciter=0; exciter<n_exciters; ++exciter)
    {
      for (unsigned int sensor=0; sensor<n_sensors; ++sensor)
      {
        Tensor<1,dim> h0_re;
        Tensor<1,dim> h0_im;
        for (unsigned int d=0;d<dim; ++d)
        {
          h0_re[d]=all_h0_at_centres[exciter][sensor](d);
          h0_im[d]=all_h0_at_centres[exciter][sensor](d+dim);
        }
        measured_voltages_re[exciter](sensor) = all_sensor_directions[sensor]*
        (h0_re + integral_factor*all_hdelta_re[exciter][sensor]);
        measured_voltages_im[exciter](sensor) = all_sensor_directions[sensor]*
        (h0_im + integral_factor*all_hdelta_im[exciter][sensor]);

      }
    }
  }
  template void alternateVoltage<3, DoFHandler<3>>;
  
} // END namespace Voltage