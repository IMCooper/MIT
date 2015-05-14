#include <inversesolver.h>

using namespace dealii;

namespace InverseSolver
{
  template <int dim>
  void process_mesh(Triangulation<dim> &triangulation)
  {
    // This routine processes a mesh so that it is ready for use
    // with the InverseSolver class.
    //
    // This requires renumber the cell material_ids. 
    // It also calculates the number of voxels and coil combinations
    // and updates the data in InverseProblemData.
    
    // The original mesh must contain a recovery and non-recovery region typically
    // this will be a cube inside a cylinder, with the cube being for recovery.
    // The default non-recovery materialID is 0, or it may be set in the input file (TODO)
    // The default recovery region is 1, or it may be set in the input file (TODO)
    
    // Loop over all cells and set the material_id of
    // any non-recovery region to 0, and
    // set each recovery region cell incrementally.
    typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
    
    unsigned int recovery_cell_count=0;;
    for (; cell!=endc; ++cell)
    {
      if (cell->material_id() == InverseProblemData::recovery_region_id)
      {
        ++recovery_cell_count;
        cell->set_material_id(recovery_cell_count);
      }   
      else if (cell->material_id() == InverseProblemData::background_region_id)
      {
        cell->set_material_id(0);
      }
      else
      {
        // TODO: use deal error properly.
        std::cout << "ERROR: Input mesh not compatibile with InverseSolver" << std::endl; 
        std::cout << "Expected material_id: "
        << (unsigned int)InverseProblemData::recovery_region_id
        << " or "
        << (unsigned int)InverseProblemData::background_region_id << std::endl;
        std::cout << "Found material_id: " << (unsigned int)cell->material_id() << std::endl;
        
      }
    }  
    
    // number of voxels is given by
    InverseProblemData::n_voxels = recovery_cell_count;
    
    // Now refine the outer mesh
    cell = triangulation.begin_active();
    for (; cell!=endc; ++cell)
    {
//       if (cell->material_id() == 0)
//       {
        cell->set_refine_flag();
//       }                         
    }
    triangulation.execute_coarsening_and_refinement ();    
  }
  
//   template <int dim>
  void initialise_material_parameters()
  {
    // Once the number of voxels is known, we can then initialise
    // the recovery parameters. For now we begin with the background value
    // for every voxel.
    // Note the number of materials is InverseProblemData::n_voxels+1
    
    // first resize the vectors:
    EquationData::param_mur.reinit(InverseProblemData::n_voxels+1);
    EquationData::param_sigma.reinit(InverseProblemData::n_voxels+1);
    EquationData::param_epsilon.reinit(InverseProblemData::n_voxels+1);
    EquationData::param_kappa_re.reinit(InverseProblemData::n_voxels+1);
    EquationData::param_kappa_im.reinit(InverseProblemData::n_voxels+1);
    
    InverseProblemData::initial_sigma.reinit(InverseProblemData::n_voxels+1);
    
    // Now fill with the background values.
    // Note:
    // kappa = -omega^2*epr + i*omega*sigma;
    // i.e. kappa_re = -omega^2
    //      kappa_im = omega*sigma
    
    // Set the outer region to be vacuum values:
    EquationData::param_mur(0) = 1.0;
    EquationData::param_sigma(0) = 0.0;
    EquationData::param_epsilon(0) = 0.0; // zero as using eddy current approx.
    for (unsigned int i=1; i<InverseProblemData::n_voxels+1; ++i)
    {
      EquationData::param_mur(i) = EquationData::param_mur_background;
      EquationData::param_sigma(i) = EquationData::param_sigma_background;
      EquationData::param_epsilon(i) = EquationData::param_epsilon_background;
      
      // TESTING:
      // Add conducting region 
      /*
      if ( (i > 11 + 25 && i < 15 + 25) || (i > 11 + 50 && i < 15 + 50) || (i > 11 + 75 && i < 15 + 75) )
      {
        EquationData::param_sigma(i) = EquationData::param_sigma_conducting;
      }
      */
    }
    // Now fill kappa:
    for (unsigned int i=0; i<InverseProblemData::n_voxels+1; ++i)
    {
      if (EquationData::param_sigma(i) > EquationData::param_regularisation)
      {
        EquationData::param_kappa_re(i) = 0.0;//EquationData::param_omega*EquationData::param_sigma(i)*EquationData::constant_mu0;
        EquationData::param_kappa_im(i) = EquationData::param_sigma(i)*EquationData::param_omega*EquationData::constant_mu0;
      }
      else
      {
        EquationData::param_kappa_re(i) = EquationData::param_regularisation*EquationData::param_omega*EquationData::constant_mu0;
        EquationData::param_kappa_im(i) = 0.0;//EquationData::param_regularisation*EquationData::constant_mu0*EquationData::param_omega;
      }
    }
    InverseProblemData::initial_sigma = EquationData::param_sigma;
  }
  
//   template <int dim>
  void setup_coils()
  {
    // Uses the information in ExcitationCoilData and SensorCoilData
    // to set up the coil positions and directions for each (and fill the 
    // storage within each namespace).
    
    //TODO: could later add in a tilt/rotation angle.. for now this isn't needed.
    
    // Excitation coils:
    // resize the vectors:
    ExcitationCoilData::coil_position.resize(ExcitationCoilData::number_of_coils);
    ExcitationCoilData::coil_direction.resize(ExcitationCoilData::number_of_coils);   
    // Calculate the positions/directions:
    double excitation_direction_factor = numbers::PI*ExcitationCoilData::coil_radius*ExcitationCoilData::coil_radius;
    for (unsigned int coil_i=0; coil_i < ExcitationCoilData::number_of_coils; ++coil_i)
    { 
      double angle = ExcitationCoilData::array_angle + double(coil_i)*2.0*numbers::PI/double(ExcitationCoilData::number_of_coils);
      Point<3> temp_point (ExcitationCoilData::array_radius*cos(angle) + ExcitationCoilData::array_centre[0],
                             ExcitationCoilData::array_radius*sin(angle) + ExcitationCoilData::array_centre[1],
                             ExcitationCoilData::array_centre[2]);
      
      Tensor<1,3> temp_direction({excitation_direction_factor*cos(angle),
                                    excitation_direction_factor*sin(angle),
                                    0.0});
      
      ExcitationCoilData::coil_position[coil_i] = temp_point;
      ExcitationCoilData::coil_direction[coil_i] = temp_direction;
    }
    
    // Sensor coils:
    // resize the vectors:
    SensorCoilData::coil_position.resize(SensorCoilData::number_of_coils);
    SensorCoilData::coil_direction.resize(SensorCoilData::number_of_coils);
    // Calculate the positions/directions:
    double sensor_direction_factor = numbers::PI*SensorCoilData::coil_radius*SensorCoilData::coil_radius;
    for (unsigned int coil_i=0; coil_i < SensorCoilData::number_of_coils; ++coil_i)
    {
      double angle = SensorCoilData::array_angle + double(coil_i)*2.0*numbers::PI/double(SensorCoilData::number_of_coils);
      Point<3> temp_point (SensorCoilData::array_radius*cos(angle) + SensorCoilData::array_centre[0],
                             SensorCoilData::array_radius*sin(angle) + SensorCoilData::array_centre[1],
                             SensorCoilData::array_centre[2]);
      
      Tensor<1,3> temp_direction({excitation_direction_factor*cos(angle),
                                  excitation_direction_factor*sin(angle),
                                  0.0});
      
      SensorCoilData::coil_position[coil_i] = temp_point;
      SensorCoilData::coil_direction[coil_i] = temp_direction;
    }
  }
  
  template <int dim>
  void initialise_inverse_problem(Triangulation<dim> &triangulation)
  {
    // MUST CALL THIS BEFORE THE INVERSESOLVER CLASS
    // This routine will setup the inverse problem, calling
    // the routines above to setup the mesh, the coil locations/directions
    // and will initialise all of the required variables/data.    
    process_mesh(triangulation);
    initialise_material_parameters();
    setup_coils();
    
    // Calculate the coil combinations - used for the sensitivity matrix size.
    // the number of coil combinations is
    // (the number of excitation coils) * (the number of sensor coils)    
    InverseProblemData::n_coil_combinations
      = SensorCoilData::number_of_coils*ExcitationCoilData::number_of_coils;
    // Uncomment to avoid coil_i = coil_j combinations (& add the use_all part to input sections - TODO)
//     if (InverseProblemData::use_all)
//     {
//       InverseProblemData::n_coil_combinations
//       = SensorCoilData::number_of_coils*ExcitationCoilData::number_of_coils;
//     }
//     else
//     {
//       InverseProblemData::n_coil_combinations
//       = ExcitationCoilData::number_of_coils*(SensorCoilData::number_of_coils-1);
//     }
  }
  template void initialise_inverse_problem<3>;
  
  // InverseSolver class functions:  
  template<int dim, class DH>
  InverseSolver<dim, DH>::InverseSolver( std::vector<Vector<double>> &solutions_in,
                                         FiniteElement<dim> &fe)
  :
  // Constructor: attach dof handler and copy solution vectors.
  solutions(solutions_in),
  fe(&fe)
  {
  }
  
  template<int dim, class DH>
  void InverseSolver<dim, DH>::assemble_sensitivity_matrix(const DH &dof_handler)
  {
    // Routine to assemble the sensitivity matrix for the MIT problem
    // This will be stored in the private FullMatrix, sensitivity_matrix.
    
    // This routine calculates the sensitivity matrix for the given
    // set of solutions in the constructor and uses the DoFHandler to 
    // calculate the entries via the finite element
    //
    // The entries are given by:
    // J_ij = \int_vox_{i}} E_{a} \cdot E_{b} dV
    // 
    // a = the exciter coil index,
    // b = the sensor coil index,
    // i = b + a*number_of_sensor_coils;
    
    // resize matrices:
    sensitivity_matrix_re.reinit(InverseProblemData::n_coil_combinations,
                                 InverseProblemData::n_voxels);
    sensitivity_matrix_im.reinit(InverseProblemData::n_coil_combinations,
                                 InverseProblemData::n_voxels);
    
    QGauss<dim>  quadrature_formula((fe->degree+1)*2+1);
    const unsigned int n_q_points = quadrature_formula.size();
    
    const unsigned int dofs_per_cell = fe->dofs_per_cell;
    
    FEValues<dim> fe_values (*fe, quadrature_formula,
                             update_values |
                             update_quadrature_points  |  update_JxW_values);
    
    // Extractors to real and imaginary parts
      const FEValuesExtractors::Vector E_re(0);
      const FEValuesExtractors::Vector E_im(dim);
    
    // Storage for computed solution on a cell at each quadrature point
    std::vector<std::vector<Tensor<1,dim>>> local_solutions_re(ExcitationCoilData::number_of_coils,
                                                               std::vector<Tensor<1,dim>>(n_q_points));
    std::vector<std::vector<Tensor<1,dim>>> local_solutions_im(ExcitationCoilData::number_of_coils,
                                                               std::vector<Tensor<1,dim>>(n_q_points));
    
    sensitivity_matrix_re=0;
    sensitivity_matrix_im=0;
    
    typename DH::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    
    // Loop over all cells
    for (; cell!=endc; ++cell)
    {
      // take note of which voxel we are in (via material_id)
      // If material_id is not in the recovery region, then skip.
      // The background region is hard-coded to be 0 in the setup region
      if (cell->material_id() > 0)
      {
        // Since 0 is the background region, take 1 from the material_id to get the voxel index
        // i.e. material_id 1 is voxel 0 & so on.
        unsigned int voxel_index = (unsigned int)cell->material_id() - 1;
        
        fe_values.reinit(cell);
        
        // within this cell, calculate solution values at all quadrature points:
        for (unsigned int exciter=0; exciter < ExcitationCoilData::number_of_coils; ++exciter)
        {
          fe_values[E_re].get_function_values(solutions[exciter], local_solutions_re[exciter]);
          fe_values[E_im].get_function_values(solutions[exciter], local_solutions_im[exciter]);
        }
        
        // begin quadrature loop:      
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
          for (unsigned int exciter=0; exciter<ExcitationCoilData::number_of_coils; ++exciter)
          {
            for (unsigned int sensor=0; sensor<SensorCoilData::number_of_coils; ++ sensor)
            {
              unsigned int i_index = sensor + exciter*(SensorCoilData::number_of_coils);
              
              sensitivity_matrix_re(i_index, voxel_index)
              += ( local_solutions_re[exciter][q_point]*local_solutions_re[sensor][q_point]
              - local_solutions_im[exciter][q_point]*local_solutions_im[sensor][q_point] )*fe_values.JxW(q_point);
              
              sensitivity_matrix_im(i_index, voxel_index)
              += ( local_solutions_re[exciter][q_point]*local_solutions_im[sensor][q_point]
              + local_solutions_im[exciter][q_point]*local_solutions_re[sensor][q_point] )*fe_values.JxW(q_point);            
            }
          }
        }
      }
    }
    // Due to the boundary conditions, we are solving the "A-based" formulation, so we must
    // multiply the solution by -i*omega to get the electric field.
    // Create a factor to multiply the sensitivity matrix by to get this:
    double sensitivity_factor = -EquationData::param_omega*EquationData::param_omega;
    sensitivity_matrix_re *= sensitivity_factor;
    sensitivity_matrix_im *= sensitivity_factor;

    // TESTING: print out sens matrices (re/im).
    {
    std::stringstream matrix_filename;
    matrix_filename << "sens_matrix_re_it"
    << InverseProblemData::iteration_count
    << "_p" << fe->degree-1 << ".txt";
    std::ofstream matrix_file(matrix_filename.str());
    matrix_file.precision(32);
    for (unsigned int i=0; i<InverseProblemData::n_coil_combinations; ++i)
    {
      
      for (unsigned int j=0; j<InverseProblemData::n_voxels; ++j)
      {
        matrix_file << sensitivity_matrix_re(i,j) << " ";
      }
      matrix_file << std::endl;
    }
    matrix_file.close();
    }
    {
    std::stringstream matrix_filename;
    matrix_filename << "sens_matrix_im_it"
    << InverseProblemData::iteration_count
    << "_p" << fe->degree-1 << ".txt";
    std::ofstream matrix_file(matrix_filename.str());
    matrix_file.precision(32);
    for (unsigned int i=0; i<InverseProblemData::n_coil_combinations; ++i)
    {
      
      for (unsigned int j=0; j<InverseProblemData::n_voxels; ++j)
      {
        matrix_file << sensitivity_matrix_im(i,j) << " ";
      }
      matrix_file << std::endl;
    }
    matrix_file.close();
    }
  }
  
  template<int dim, class DH>
  void InverseSolver<dim, DH>::assemble_sensitivity_rhs(const DH &dof_handler)
  {
    // Construct the RHS of the sensitivity system. This is:
    // -1/2*(J^{H}*V + J^{T}*V^{*}).
    //
    // where J is the sensitivity matrix calculated by assemble_sensitivity_matrix()
    // and V is the voltage difference between the image measured voltages and the
    // current simulated voltages.
    //
    // ^{*} denotes complex conjugate and ^{H} denotes hermitian (conjugate transpose).
    //
    // Care must be taken to match the exciter/sensor coil combinations of the voltage
    // to that which was used in the assemble_sensitivity_matrix() routine.    
    
    // resize the storage for sens rhs:
    sensitivity_rhs_re.reinit(InverseProblemData::n_voxels);
    sensitivity_rhs_im.reinit(InverseProblemData::n_voxels);
    
    // First calculate the simulated voltages
    // TODO: Consider moving the voltage calculation outside of the inverse solver,
    //       and instead passing the voltage differences to this function ???
    std::vector<Vector<double>>
    simulated_voltages_re(ExcitationCoilData::number_of_coils,
                          Vector<double> (SensorCoilData::number_of_coils));
    std::vector<Vector<double>>
    simulated_voltages_im(ExcitationCoilData::number_of_coils,
                          Vector<double> (SensorCoilData::number_of_coils));
    for (unsigned int exciter=0; exciter<ExcitationCoilData::number_of_coils; ++exciter)
    {
      for (unsigned int sensor=0; sensor<SensorCoilData::number_of_coils; ++sensor)
      {
        Voltage::Voltage<dim> measure_voltage(SensorCoilData::coil_position[sensor],
                                              SensorCoilData::coil_direction[sensor]);
      
        measure_voltage.calculateVoltage(solutions[exciter],
                                         dof_handler,
                                         simulated_voltages_re[exciter](sensor),
                                         simulated_voltages_im[exciter](sensor));
      }
    }
    
    // calc number of exciter/sensor combinations
    // (i.e. we have N exciters each with M sensors = N*M.
    Vector<double> voltage_differences_re(InverseProblemData::n_coil_combinations);
    Vector<double> voltage_differences_im(InverseProblemData::n_coil_combinations);
    
    std::vector<Vector<double>> temp_voltage_diffs_re(ExcitationCoilData::number_of_coils,
                                                      Vector<double> (SensorCoilData::number_of_coils));
    std::vector<Vector<double>> temp_voltage_diffs_im(ExcitationCoilData::number_of_coils,
                                                      Vector<double> (SensorCoilData::number_of_coils));
    for (unsigned int exciter=0; exciter<ExcitationCoilData::number_of_coils; ++exciter)
    {
      for (unsigned int sensor=0; sensor<SensorCoilData::number_of_coils; ++sensor)
      {
        unsigned int i_index = sensor + exciter*(SensorCoilData::number_of_coils);
        
        voltage_differences_re[i_index]
        = simulated_voltages_re[exciter](sensor) - InverseProblemData::measured_voltages_re[exciter](sensor);

        voltage_differences_im[i_index]
        = simulated_voltages_im[exciter](sensor) - InverseProblemData::measured_voltages_im[exciter](sensor);
        
        //TESTING:
//         if (abs(voltage_differences_im[i_index]) > 0 || abs(voltage_differences_re[i_index]) > 0)
//         {
//           std::cout << exciter << " " << sensor << " " << i_index << " "
//           << voltage_differences_re[i_index] << " " << voltage_differences_im[i_index] << std::endl;
//         std::cout << simulated_voltages_re[exciter](sensor) << " " << simulated_voltages_im[exciter](sensor) << std::endl;
//         std::cout << InverseProblemData::measured_voltages_re[exciter](sensor) << " " << InverseProblemData::measured_voltages_im[exciter](sensor) << std::endl;
//         }
        
        
        temp_voltage_diffs_re[exciter](sensor) = voltage_differences_re[i_index];
        temp_voltage_diffs_im[exciter](sensor) = voltage_differences_im[i_index];
      }
    }
    {
      // TESTING:
      {
      std::stringstream filename;
      filename << "voltage_diffs"
      << "_p" << dof_handler.get_fe().degree-1  << ".txt";
      OutputTools::output_voltages(temp_voltage_diffs_re,
                                   temp_voltage_diffs_im,
                                   filename.str());      
      }
      {
      std::stringstream filename;      
      filename << "simulated_volts"
      << "_p" << dof_handler.get_fe().degree-1  << ".txt";
      OutputTools::output_voltages(simulated_voltages_re,
                                   simulated_voltages_im,
                                   filename.str());
      }
      {
      std::stringstream filename;
      filename <<"measured_volts"
      << "_p" << dof_handler.get_fe().degree-1  << ".txt";
      OutputTools::output_voltages(InverseProblemData::measured_voltages_re,
                                   InverseProblemData::measured_voltages_im,
                                   filename.str());    
      }
    }
    // TODO:: store the voltage_differences as they make up the residual of the inverse iteration.
    
    // calc RHS:
    // Note, complex so, J = sens matrix, V = voltage differences:
    //
    // ^{H} = hermitian. ^{*} = complex conjugate.
    //
    // J^{H}*V = (J_re^{T} - i*J_im^{T})*(V_re + i*V_im)
    //         = (J_re^{T}*V_re + J_im^{T}*V_im) + i*(J_re^{T}*V_im - J_im^{T}*V_re).
    //
    // J^{T}*V^{*} = (J_re^{T} + i*J_im^{T})*(V_re - i*V_im)
    //             = (J_re^{T}*V_re + J_im^{T}*V_im) + i*(J_im^{T}*V_re - J_re^{T}*V_im).
    //
    // So RHS is..
    // real part: 1/2*( (J_re^{T}*V_re + J_im^{T}*V_im) + (J_re^{T}*V_re + J_im^{T}*V_im) )
    //            = J_re^{T}*V_re + J_im^{T}*V_im
    //
    // imag part: 1/2*(J_re^{T}*V_im - J_im^{T}*V_re) + (J_im^{T}*V_re - J_re^{T}*V_im) = 0.
    
    // Real part:
    sensitivity_matrix_re.Tvmult(sensitivity_rhs_re,
                                 voltage_differences_re); // copies over sensitivity_rhs_re
    sensitivity_matrix_im.Tvmult_add(sensitivity_rhs_re,
                                     voltage_differences_im); // adds to sensitivity_rhs_re
    
    // TESTING: output residual:
    double temp_residual=0;
    for (unsigned int i=0; i<InverseProblemData::n_coil_combinations; ++i)
    {
      temp_residual += voltage_differences_re[i]*voltage_differences_re[i] + voltage_differences_im[i]*voltage_differences_im[i];
    }
      
      std::cout << "inverse residual: " << sqrt(temp_residual) << " " << EquationData::param_sigma(1) << std::endl;

  }
  
  template<int dim, class DH>
  void InverseSolver<dim, DH>::assemble_regularisation(const DH &dof_handler)
  {
    // TODO: add regularisation matrix. For now, using identity.
    
  }
  
  template<int dim, class DH>
  void InverseSolver<dim, DH>::gauss_newton_solve(const Vector<double> &last_solution,
                                                  Vector<double> &solution_out)
  {
    // TODO: Add output of the residual (voltage differences).
    // Could include the regularisation part of the residual as well (1/2*lambda^"norm(reg_matrix*dX)) etc.
    
    // Input is the last computed solution (deltaX, not X), which is needed for regularisation
    // and the output is the computed solution.
    //
    // Note that both vectors should be of length 2*InverseProblemData::n_voxels (real & imaginary components)
    //
    //
    // Gauss-Newton iterative inverse solver with regularisation (TODO)
    // Assumes we've already called:
    // - assemble_sensitivity_matrix
    // - assemble_sensitivity_rhs
    // - assemble_regularisation (TODO)
    //
    // Update term is:
    // dX = (Re(J^{H}*J) + alpha^2*R^{T}*R)^{-1}*(sensitivity_rhs + alpha^2*R^{T}*R*X)
    //
    // where X is the parameter distribution,
    //       R is the regularisation matrix,
    //       alpha is the regularisation parameter.
    //       J is the sensitivity matrix.
    // where 
    // Note for some matrix J = J_re + i*J_im:
    // say A = J^{H}*J = (J_re^{T} - i*J_im^{T})*(J_re + i*J_im)
    //                 = (J_re^{T}*J_re + J_im^{T}*J_im) + i*(J_im^{T}*J_re - J_re^{T}*J_im)
    //
    // i.e. A_re = J_re^{T}*J_re + J_im^{T}*J_im
    //      A_im = J_im^{T}*J_re - J_re^{T}*J_im
    //
    // However, here the matrix is Re(J^{H}*J), so we only need J_re^{T}*J_re + J_im^{T}*J_im.
    
    //increment iteration count:
    ++InverseProblemData::iteration_count;
    
    
    unsigned int n_voxels = InverseProblemData::n_voxels;
    FullMatrix<double> GN_matrix(n_voxels, n_voxels);
    // Options for computing the inverse are:
    // Use GMRES as the real system arising from the complex matrix is non-symmetric.
    // Alternative is to use a direct solver via FullMatrix::invert() or gauss_jordan()
    // TODO: Use gauss jordan for now, but should experiment with performance.
    
    // Construct the GN matrix:
    // Real part only:
    sensitivity_matrix_re.Tmmult(GN_matrix,
                                 sensitivity_matrix_re); // copies over GN_matrix.
    sensitivity_matrix_im.Tmmult(GN_matrix,
                                 sensitivity_matrix_im,
                                 true); // adds result to GN_matrix.

    // For now we're using the identity for the regularisation matrix
    // so just add the regularisation parameter to the diagonal entries.
    // TODO: update to include a general regularisation matrix. For now just use identity.
    for (unsigned int i=0; i<n_voxels; ++i)
    {
      GN_matrix(i,i) += regularisation_parameter*regularisation_parameter;
    }
    
    // TESTING:
    std::stringstream matrix_filename;
    matrix_filename << "gn_matrix_it" << InverseProblemData::iteration_count << ".txt";
    std::ofstream matrix_file(matrix_filename.str());
    matrix_file.precision(32);
    for (unsigned int i=0; i<n_voxels; ++i)
    {
      
      for (unsigned int j=0; j<n_voxels; ++j)
      {
        matrix_file << GN_matrix(i,j) << " ";
      }
      matrix_file << std::endl;
    }
    matrix_file.close();
    
    // calculate the inverse: GN_matrix will now hold it.
    GN_matrix.gauss_jordan(); 
    
    // Now form the RHS:
    // For now only using identity regularisation matrix
    // TODO: Update for a general regularisation matrix.
    Vector<double> GN_rhs(n_voxels);
    for (unsigned int i=0; i<n_voxels; ++i)
    {
      // removed reg term: not sure if it should be there. Without this, we get a single step tikhonov regularisation ??
      GN_rhs(i) = sensitivity_rhs_re(i)
      + regularisation_parameter*regularisation_parameter*(EquationData::param_sigma(i+1) - InverseProblemData::initial_sigma(i+1)); 
    }
    
    //TESTING
    std::stringstream rhs_filename;
    rhs_filename << "gn_rhs_it" << InverseProblemData::iteration_count << ".txt";
    std::ofstream rhs_file(rhs_filename.str());
    rhs_file.precision(32);
    for (unsigned int i=0; i<n_voxels; ++i)
    {
      rhs_file << GN_rhs(i) << std::endl;
    }
    rhs_file.close();
    
    // resize solution:
    solution_out.reinit(n_voxels);
    // Multiply RHS by the inverted matrix:    
    GN_matrix.vmult(solution_out, GN_rhs);
  }
  template class InverseSolver<3, DoFHandler<3>>;
}