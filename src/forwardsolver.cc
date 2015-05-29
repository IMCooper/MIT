#include <forwardsolver.h>

using namespace dealii;

namespace ForwardSolver
{
  // CLASS EDDY CURRENT MEMBERS:
  template <int dim, class DH>
  EddyCurrent<dim, DH>::EddyCurrent(DH &dof_handler,
                                    const FiniteElement<dim> &fe,
                                    const bool direct_flag)
  :
  fe(&fe),
  mapping(&StaticMappingQ1<dim>::mapping),
  direct(direct_flag)
  {
    
    // Constructor for the class:
    // 
    // Input:
    // - DoFHandler with a triangulation and finite element attached (can used dof_handler.initialise).
    // - Finite element, which should be a FESystem containing two blocks/copies of an FE_Nedelec element.
    // 
    // Output:
    // - The modified DoFHandler.
    
    constructor_setup(dof_handler,
                      direct_flag);
  }
  template <int dim, class DH>
  EddyCurrent<dim, DH>::EddyCurrent(const Mapping<dim> &mapping_in,
                                    DH &dof_handler,
                                    const FiniteElement<dim> &fe,
                                    const bool direct_flag)
  :
  fe(&fe),
  mapping(&mapping_in),
  direct(direct_flag)
  {
    constructor_setup(dof_handler,
                      direct_flag);
  }
  
  template<int dim, class DH>
  void EddyCurrent<dim, DH>::constructor_setup(DH &dof_handler,
                                          const bool direct_flag)
                                          
  {
    // Code common to either version of the contructor
    // 
    // This will distribute the degrees of freedom for the input FiniteElement and then re-order
    // the DoF numbering into blocks of lower and higher order basis functions.
    //
    // Finally, the constraints and sparsity pattern for the matrix is created using a dummy boundary function.
    // If the input DoFHandler was attached to a triangulation with non-Dirichlet BCs (boundary_id > 0) then
    // the constraints will only contain hanging node constraints (if present).
    
    // set private direct solver flag:
    direct = direct_flag;

    p_order = fe->degree - 1;
    
    //TODO: work out what the order should be.
    // Some quick tests suggest it should be 2*N+4..
    // this is when the voltages didn't change against the next order up (2*N+5).
    quad_order = 2*(fe->degree) + 1;
    
    
    // Setup for FEM:
    system_matrix.clear();
    
    // Renumber by lower/higher order
    std::vector<unsigned int> reorder_counts;
    MyDoFRenumbering::by_dimension<dim, DoFHandler<dim>>(dof_handler, reorder_counts);
    
    n_lowest_order_dofs = reorder_counts[0];
    n_higher_order_edge_gradients_dofs = reorder_counts[1];
    n_higher_order_face_gradients_dofs = reorder_counts[2];
    n_higher_order_cell_gradients_dofs = reorder_counts[3];
    n_higher_order_face_nongradients_dofs = reorder_counts[4];
    n_higher_order_cell_nongradients_dofs = reorder_counts[5];
    
    n_higher_order_dofs
    = n_higher_order_edge_gradients_dofs
    + n_higher_order_face_gradients_dofs
    + n_higher_order_face_nongradients_dofs
    + n_higher_order_cell_gradients_dofs
    + n_higher_order_cell_nongradients_dofs;

    n_higher_order_gradient_dofs
    = n_higher_order_edge_gradients_dofs
    + n_higher_order_face_gradients_dofs
    + n_higher_order_cell_gradients_dofs;
    
    n_higher_order_non_gradient_dofs
    = n_higher_order_face_nongradients_dofs
    + n_higher_order_cell_nongradients_dofs;
    
    int remaining_dofs = dof_handler.n_dofs()
    - n_lowest_order_dofs
    - n_higher_order_edge_gradients_dofs
    - n_higher_order_face_gradients_dofs
    - n_higher_order_face_nongradients_dofs
    - n_higher_order_cell_gradients_dofs
    - n_higher_order_cell_nongradients_dofs;
    
    if (remaining_dofs !=0)
    {
      std::cout << std::endl << "WARNING! Renumbering did not find all DoFs! " << remaining_dofs << " were not renumbered." << std::endl;
    }
    
    // Setup initial constraints in order to pre-construct the sparsity pattern:
    constraints.clear ();
    
    // Hanging nodes come second according to project_boundary_values_curl_conforming_l2 documentation:
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);
    
    // FE_Nedelec boundary constraints: note we just a zero function as we have not passed the actual function yet.
    // TODO: verify that this does not change anything !!
    // Real part (begins at 0):
    VectorTools::project_boundary_values_curl_conforming_l2 (dof_handler,
                                                             0,
                                                             ZeroFunction<dim> (dim+dim),
                                                             0,
                                                             constraints,
                                                             *mapping);
    // Imaginary part (begins at dim):
    VectorTools::project_boundary_values_curl_conforming_l2 (dof_handler,
                                                             dim,
                                                             ZeroFunction<dim> (dim+dim),
                                                             0,
                                                             constraints,
                                                             *mapping);

    // Add constraints on gradient-based DoFs outside of the conducting region,
    // i.e. in cells with material_id 0.
    // These DoFs are automatically zero in this region because the curl-curl operator
    // applied to a gradient-based DoF is zero and the coefficient of the mass operator
    // is effectively zero.
    // Note: we must avoid constraining DoFs on the suface of the conductor, so first
    // we make a list of DoFs inside the conductor
    if (p_order > 0)
    {
      std::vector<bool> dof_in_conductor (dof_handler.n_dofs());
      std::vector<types::global_dof_index> local_dof_indices (fe->dofs_per_cell);
      
      typename DoFHandler<dim>::active_cell_iterator cell, endc;
      endc = dof_handler.end();
      cell = dof_handler.begin_active();
      
      for (;cell!=endc; ++cell)
      {
        // TODO: update with the conductor material id from stored data.
        if (cell->material_id() == 0)
        {
          cell->get_dof_indices (local_dof_indices);
          for (unsigned int i = 0; i<fe->dofs_per_cell; ++i)
          {
            dof_in_conductor[local_dof_indices[i]] = true;
          }
        }
      }
      // Now loop through all cells in the non-conducting region and set gradient DoFs to be zero,
      // on the condition that they are not flagged by dof_in_conductor.
      cell = dof_handler.begin_active();
      for (;cell!=endc; ++cell)
      {
        // TODO: update with the conductor material id from stored data.
        if (cell->material_id() == 1)
        {
          cell->get_dof_indices (local_dof_indices);
          for (unsigned int i = 0; i<fe->dofs_per_cell; ++i)
          {
            if (local_dof_indices[i] >= n_lowest_order_dofs
              && local_dof_indices[i] < n_higher_order_gradient_dofs
              && !dof_in_conductor[local_dof_indices[i]])
            {
              constraints.add_line(local_dof_indices[i]);
            }
          }
        }
      }
    }
    
    constraints.close ();

    if (p_order == 0)
    {
      BlockDynamicSparsityPattern csp (1,1);
      csp.block(0,0).reinit (n_lowest_order_dofs, n_lowest_order_dofs);
      csp.collect_sizes();
      
      DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);
      
      sparsity_pattern.copy_from(csp);
    }
    else
    {
      // 2 by 2 block:
//       BlockDynamicSparsityPattern csp (2,2);
//       
//       csp.block(0,0).reinit (n_lowest_order_dofs, n_lowest_order_dofs);
//       csp.block(0,1).reinit (n_lowest_order_dofs, n_higher_order_dofs);
//       csp.block(1,0).reinit (n_higher_order_dofs, n_lowest_order_dofs);
//       csp.block(1,1).reinit (n_higher_order_dofs, n_higher_order_dofs);
//       csp.collect_sizes();
      
      // 3 by 3 block
      // Block1: Lowest order
      // Block2: gradient-based higher order
      // Block3: non-gradient-based higher order
      BlockDynamicSparsityPattern csp(3,3);
      csp.block(0,0).reinit (n_lowest_order_dofs,n_lowest_order_dofs);
      csp.block(0,1).reinit (n_lowest_order_dofs, n_higher_order_gradient_dofs);
      csp.block(0,2).reinit (n_lowest_order_dofs, n_higher_order_non_gradient_dofs);
      
      csp.block(1,0).reinit (n_higher_order_gradient_dofs, n_lowest_order_dofs);
      csp.block(1,1).reinit (n_higher_order_gradient_dofs, n_higher_order_gradient_dofs);
      csp.block(1,2).reinit (n_higher_order_gradient_dofs, n_higher_order_non_gradient_dofs);
      
      csp.block(2,0).reinit (n_higher_order_non_gradient_dofs, n_lowest_order_dofs);
      csp.block(2,1).reinit (n_higher_order_non_gradient_dofs, n_higher_order_gradient_dofs);
      csp.block(2,2).reinit (n_higher_order_non_gradient_dofs, n_higher_order_non_gradient_dofs);
      csp.collect_sizes();
      //Full 4 by 4 block:
      /*
       *   BlockDynamicSparsityPattern csp (4,4);
       *   
       *   csp.block(0,0).reinit (n_lowest_order_dofs, n_lowest_order_dofs);
       *   csp.block(0,1).reinit (n_lowest_order_dofs, n_higher_order_edge_dofs);
       *   csp.block(0,2).reinit (n_lowest_order_dofs, n_higher_order_face_dofs);
       *   csp.block(0,3).reinit (n_lowest_order_dofs, n_higher_order_cell_dofs);
       *   
       *   csp.block(1,0).reinit (n_higher_order_edge_dofs, n_lowest_order_dofs);
       *   csp.block(1,1).reinit (n_higher_order_edge_dofs, n_higher_order_edge_dofs);
       *   csp.block(1,2).reinit (n_higher_order_edge_dofs, n_higher_order_face_dofs);
       *   csp.block(1,3).reinit (n_higher_order_edge_dofs, n_higher_order_cell_dofs);
       *   
       *   csp.block(2,0).reinit (n_higher_order_face_dofs, n_lowest_order_dofs);
       *   csp.block(2,1).reinit (n_higher_order_face_dofs, n_higher_order_edge_dofs);
       *   csp.block(2,2).reinit (n_higher_order_face_dofs, n_higher_order_face_dofs);
       *   csp.block(2,3).reinit (n_higher_order_face_dofs, n_higher_order_cell_dofs);
       *   
       *   csp.block(3,0).reinit (n_higher_order_cell_dofs, n_lowest_order_dofs);
       *   csp.block(3,1).reinit (n_higher_order_cell_dofs, n_higher_order_edge_dofs);
       *   csp.block(3,2).reinit (n_higher_order_cell_dofs, n_higher_order_face_dofs);
       *   csp.block(3,3).reinit (n_higher_order_cell_dofs, n_higher_order_cell_dofs);
       *   csp.collect_sizes();
       */
      
      DoFTools::make_sparsity_pattern(dof_handler,  csp, constraints, false);
      
      sparsity_pattern.copy_from(csp);
    }
    system_matrix.reinit (sparsity_pattern);
    if (!direct)
    {
      system_preconditioner.reinit (sparsity_pattern);
    }
    
    if (p_order==0)
    {
      solution.reinit(1);
      solution.block(0).reinit (n_lowest_order_dofs);
      solution.collect_sizes();
      
      system_rhs.reinit(1);
      system_rhs.block(0).reinit(n_lowest_order_dofs);
      system_rhs.collect_sizes();
    }
    else
    {
      // 2 by 2 blocks:
//       solution.reinit(2);
//       solution.block(0).reinit (n_lowest_order_dofs);
//       solution.block(1).reinit (n_higher_order_dofs);
//       solution.collect_sizes();
//       
//       system_rhs.reinit(2);
//       system_rhs.block(0).reinit (n_lowest_order_dofs);
//       system_rhs.block(1).reinit (n_higher_order_dofs);
//       system_rhs.collect_sizes();
      
      // 3 by 3 blocks:
      solution.reinit(3);
      solution.block(0).reinit (n_lowest_order_dofs);
      solution.block(1).reinit (n_higher_order_gradient_dofs);
      solution.block(2).reinit (n_higher_order_non_gradient_dofs);
      solution.collect_sizes();
      
      system_rhs.reinit(3);
      system_rhs.block(0).reinit (n_lowest_order_dofs);
      system_rhs.block(1).reinit (n_higher_order_gradient_dofs);
      system_rhs.block(2).reinit (n_higher_order_non_gradient_dofs);
      system_rhs.collect_sizes();      
      
      // 4 by 4 blocks:
      /*
       * solution.reinit(4);
       * solution.block(0).reinit (n_lowest_order_dofs);
       * solution.block(1).reinit (n_higher_order_edge_dofs);
       * solution.block(2).reinit (n_higher_order_face_dofs);
       * solution.block(3).reinit (n_higher_order_cell_dofs);
       * solution.collect_sizes();
       * 
       * system_rhs.reinit(4);
       * system_rhs.block(0).reinit (n_lowest_order_dofs);
       * system_rhs.block(1).reinit (n_higher_order_edge_dofs);
       * system_rhs.block(2).reinit (n_higher_order_face_dofs);
       * system_rhs.block(3).reinit (n_higher_order_cell_dofs);
       * system_rhs.collect_sizes();
       */
    }
  }
  
  template <int dim, class DH>
  EddyCurrent<dim, DH>::~EddyCurrent ()
  {
    // Deconstructor: need to delete the pointers to preconditioners
    // TODO: Could this be better handled with a smart pointer?
    
    if (!direct)
    {
      delete preconditioner;
    }
    
  }
  
  template <int dim, class DH>
  void EddyCurrent<dim, DH>::compute_constraints (const DH &dof_handler,
                                                  const curlFunction<dim> &boundary_function)
  {
    // This function will compute the constraints for hanging nodes
    // and the given boundary function.
    // 
    // This function should be called every time the boundary condition has been changed.
    // 
    // The function also assumes that the DoFHandler has been correctly set up.
    // i.e. A triangulation and finite element have been attached and distribute_dofs has been called.
    
    constraints.clear ();
    
    // Hanging nodes first so that Dirichlet conditions supercede hanging node conditions.
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);
    
    // FE_Nedelec boundary condition:
    // Real part (begins at 0):
    VectorTools::project_boundary_values_curl_conforming_l2 (dof_handler,
                                                             0,
                                                             boundary_function,
                                                             0,
                                                             constraints,
                                                             *mapping);
    // Imaginary part (begins at dim):
    VectorTools::project_boundary_values_curl_conforming_l2 (dof_handler,
                                                             dim,
                                                             boundary_function,
                                                             0,
                                                             constraints,
                                                             *mapping);
    
    // Add constraints on gradient-based DoFs outside of the conducting region,
    // i.e. in cells with material_id 0.
    // These DoFs are automatically zero in this region because the curl-curl operator
    // applied to a gradient-based DoF is zero and the coefficient of the mass operator
    // is effectively zero.
    // Note: we must avoid constraining DoFs on the suface of the conductor, so first
    // we make a list of DoFs inside the conductor
    if (p_order > 0)
    {
      std::vector<bool> dof_in_conductor (dof_handler.n_dofs());
      std::vector<types::global_dof_index> local_dof_indices (fe->dofs_per_cell);
      
      typename DoFHandler<dim>::active_cell_iterator cell, endc;
      endc = dof_handler.end();
      cell = dof_handler.begin_active();
      
      for (;cell!=endc; ++cell)
      {
        // TODO: update with the conductor material id from stored data.
        if (cell->material_id() == 0)
        {
          cell->get_dof_indices (local_dof_indices);
          for (unsigned int i = 0; i<fe->dofs_per_cell; ++i)
          {
            dof_in_conductor[local_dof_indices[i]] = true;
          }
        }
      }
      // Now loop through all cells in the non-conducting region and set gradient DoFs to be zero,
      // on the condition that they are not flagged by dof_in_conductor.
      cell = dof_handler.begin_active();
      for (;cell!=endc; ++cell)
      {
        // TODO: update with the conductor material id from stored data.
        if (cell->material_id() == 1)
        {
          cell->get_dof_indices (local_dof_indices);
          for (unsigned int i = 0; i<fe->dofs_per_cell; ++i)
          {
            if (local_dof_indices[i] >= n_lowest_order_dofs
              && local_dof_indices[i] < n_higher_order_gradient_dofs
              && !dof_in_conductor[local_dof_indices[i]])
            {
              constraints.add_line(local_dof_indices[i]);
            }
          }
        }
      }
    }
    
    constraints.close ();
    
  }
  
  template <int dim, class DH>
  void EddyCurrent<dim, DH>::assemble_matrices (const DH &dof_handler)
  {
    /*
     * Function to assemble the system matrix.
     * 
     * Should really only need to be called once for a given
     * set of material parameters.
     */
    QGauss<dim>  quadrature_formula(quad_order);
    
    const unsigned int n_q_points = quadrature_formula.size();
    
    const unsigned int dofs_per_cell = fe->dofs_per_cell;
    
    FEValues<dim> fe_values (*mapping,
                             *fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points  |  update_JxW_values);
    
    
    // Extractors to real and imaginary parts
    // TODO: this would make things easier:
    const FEValuesExtractors::Vector E_re(0);
    const FEValuesExtractors::Vector E_im(dim);
    
    std::vector<FEValuesExtractors::Vector> vec(2);
    vec[0] = E_re;
    vec[1] = E_im;
    
    // Local cell storage:
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_preconditioner (dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    
    // Material parameters:
    double current_mur;
    double mur_inv;
    double current_kappa_re;
    double current_kappa_im;
    double current_kappa_magnitude;
    
    typename DH::active_cell_iterator cell, endc;
    endc = dof_handler.end();
    
    cell = dof_handler.begin_active();
    //  DEBUGGING VOLUME:
//     double volume_mat1=0;
//     double volume_mat0=0;
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);
      // FOR DEBUGGING: calculate volume of parts of the mesh:
//       for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
//       {
//         if (cell->material_id()==0)
//         {
//           volume_mat0+=fe_values.JxW(q_point);
//         }
//         if (cell->material_id()==1)
//         {
//           volume_mat1+=fe_values.JxW(q_point);
//         }
//       }
      
      current_mur = EquationData::param_mur(cell->material_id());
      mur_inv = 1.0/current_mur;
      current_kappa_re = EquationData::param_kappa_re(cell->material_id());
      current_kappa_im = EquationData::param_kappa_im(cell->material_id());
      // for preconditioner:
      current_kappa_magnitude = sqrt(current_kappa_im*current_kappa_im + current_kappa_re*current_kappa_re);
      
      
      // Store coefficients of the real/imaginary blocks:
      std::vector<std::vector<double>> kappa_matrix(2, std::vector<double> (2));
      kappa_matrix[0][0] = current_kappa_re;
      kappa_matrix[0][1] = -current_kappa_im;
      kappa_matrix[1][0] = current_kappa_im;
      kappa_matrix[1][1] = current_kappa_re;
      
      std::vector<std::vector<double>> kappa_matrix_precon(2, std::vector<double> (2));
      kappa_matrix_precon[0][0] = current_kappa_re;
      kappa_matrix_precon[0][1] = -current_kappa_im;
      kappa_matrix_precon[1][0] = current_kappa_im;
      kappa_matrix_precon[1][1] = current_kappa_re;
      
      cell_matrix = 0;
      cell_preconditioner = 0;
      
      // Options for assembly
      // TODO: decide which is quickest.
      // ORIGINAL WAY:
      // Loop over quad points:
//       for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
//       {
//         for (unsigned int i=0; i<dofs_per_cell; ++i)
//         {
//           const unsigned int block_index_i = fe->system_to_block_index(i).first;
//           // Construct local matrix:
//           for (unsigned int j=0; j<dofs_per_cell; ++j)
//           {
//             const unsigned int block_index_j = fe->system_to_block_index(j).first;
//             if (block_index_i == block_index_j)
//             {
//               if (block_index_i == 0)
//               {
//                 double curl_part = mur_inv*fe_values[E_re].curl(i,q_point)*fe_values[E_re].curl(j,q_point);
//                 double mass_part = fe_values[E_re].value(i,q_point)*fe_values[E_re].value(j,q_point);
//                 
//                 cell_matrix(i,j) += ( curl_part + current_kappa_re*mass_part )*fe_values.JxW(q_point);
//                 
//                 cell_preconditioner(i,j) += ( curl_part + current_kappa_magnitude*mass_part )*fe_values.JxW(q_point);
//               }
//               else if (block_index_i == 1)
//               {
//                 double curl_part = mur_inv*fe_values[E_im].curl(i,q_point)*fe_values[E_im].curl(j,q_point);
//                 double mass_part = fe_values[E_im].value(i,q_point)*fe_values[E_im].value(j,q_point);
//                 
//                 cell_matrix(i,j) += ( curl_part + current_kappa_re*mass_part )*fe_values.JxW(q_point);
//
//                 cell_preconditioner(i,j) += ( curl_part + current_kappa_magnitude*mass_part )*fe_values.JxW(q_point);
//               }  
//             }
//             else
//             {
//               if (block_index_i == 0) // then block_index_j == 1
//               {
//                 cell_matrix(i,j) += -current_kappa_im*fe_values[E_re].value(i,q_point)*fe_values[E_im].value(j,q_point)*fe_values.JxW(q_point);
//               }
//               else if (block_index_i == 1) // then block_index_j == 0
//               {
//                 cell_matrix(i,j) += current_kappa_im*fe_values[E_im].value(i,q_point)*fe_values[E_re].value(j,q_point)*fe_values.JxW(q_point);
//               }
//             }
//           }
//         }
//       }
      
      // NEW WAY, 1: Better use of the blocks:
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
        if (!constraints.is_constrained(local_dof_indices[i]))
        {
          const unsigned int block_index_i = fe->system_to_block_index(i).first;
          // Construct local matrix:
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            if (!constraints.is_constrained(local_dof_indices[j]))
            {
              const unsigned int block_index_j = fe->system_to_block_index(j).first;
              double mass_part = 0;
              double curl_part = 0;
              if (block_index_i == block_index_j)
              {
                for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
                {
                  curl_part
                  += fe_values[vec[block_index_i]].curl(i, q_point)
                  *fe_values[vec[block_index_j]].curl(j, q_point)
                  *fe_values.JxW(q_point);
                  
                  mass_part
                  += fe_values[vec[block_index_i]].value(i, q_point)
                  *fe_values[vec[block_index_j]].value(j, q_point)
                  *fe_values.JxW(q_point);
                }
                cell_matrix(i,j) += mur_inv*curl_part + kappa_matrix[block_index_i][block_index_j]*mass_part;
                
                cell_preconditioner(i,j) += mur_inv*curl_part + kappa_matrix_precon[block_index_i][block_index_j]*mass_part;
              }
              else // off diagonal - curl-curl operator not needed.
              {
                for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
                {
                  mass_part
                  += fe_values[vec[block_index_i]].value(i,q_point)
                  *fe_values[vec[block_index_j]].value(j,q_point)
                  *fe_values.JxW(q_point);
                }
                cell_matrix(i,j) += kappa_matrix[block_index_i][block_index_j]*mass_part;
                
                cell_preconditioner(i,j) += kappa_matrix_precon[block_index_i][block_index_j]*mass_part;
              }
            }
          }
        }
      }
      // END NEW WAY, 1:
      // NEW WAY, 2: Better use of the blocks and remove all IF statements:
//       for (unsigned int i=0; i<dofs_per_cell; ++i)
//       {
//         const unsigned int block_index_i = fe->system_to_block_index(i).first;
//         // Construct local matrix:
//         for (unsigned int j=0; j<dofs_per_cell; ++j)
//         {
//           const unsigned int block_index_j = fe->system_to_block_index(j).first;
//           double mass_part = 0;
//           double curl_part = 0;
//           for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
//           {
//             curl_part
//             += fe_values[vec[block_index_i]].curl(i,q_point)
//             *fe_values[vec[block_index_j]].curl(j,q_point)
//             *fe_values.JxW(q_point);
// 
//             mass_part
//             += fe_values[vec[block_index_i]].value(i,q_point)
//             *fe_values[vec[block_index_j]].value(j,q_point)
//             *fe_values.JxW(q_point);
//           }
//           cell_matrix(i,j) += mur_inv*curl_part + kappa_matrix[block_index_i][block_index_j]*mass_part;
//           cell_preconditioner(i,j) += ( curl_part + mass_part );
//         }
//       }
      // END NEW WAY, 2
      
      // Distribute & constraint to global matrices:
      // Note: need to apply RHS constraints using the columns of their
      // local matrix in the other routine.
      constraints.distribute_local_to_global(cell_matrix,
                                             local_dof_indices,
                                             system_matrix);
      if (!direct)
      {
        constraints.distribute_local_to_global(cell_preconditioner,
                                               local_dof_indices,
                                               system_preconditioner);
      }
    }
    // DEBUGGING:
//     std::cout << "Volume of 0 | 1 | 0+1 is: " << volume_mat0 << " | " << volume_mat1 << " | " << volume_mat0+volume_mat1 << std::endl;
  }
  
  // Assemble the rhs for a zero RHS
  // in the governing equation.
  template <int dim, class DH>
  void EddyCurrent<dim, DH>::assemble_rhs (const DH &dof_handler,
                                           const curlFunction<dim> &boundary_function)
  {
    // NOTE: This function assumes that the RHS function is zero.
    //       There is another version below which handles a non-zero RHS.
    //
    // Function to assemble the RHS for a given boundary_function,
    // which implements dirichlet and/or neumann conditions.
    // 
    // It first updates the constraints and then computes the RHS.
    //
    // Note this will completely reset the RHS stored within the class.
    
    // Zero the RHS:
    system_rhs = 0;
    
    compute_constraints(dof_handler,
                        boundary_function);
    
    QGauss<dim>  quadrature_formula(quad_order);
    const unsigned int n_q_points = quadrature_formula.size();
    
    QGauss<dim-1> face_quadrature_formula(quad_order);
    const unsigned int n_face_q_points = face_quadrature_formula.size();
    
    const unsigned int dofs_per_cell = fe->dofs_per_cell;

    // Needed to calc the local matrix for distribute_local_to_global
    // Note: only need the columns of the constrained entries.
    FEValues<dim> fe_values (*mapping,
                             *fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points  |  update_JxW_values);
    
    FEFaceValues<dim> fe_face_values(*mapping,
                                     *fe, face_quadrature_formula,
                                     update_values | update_quadrature_points |
                                     update_normal_vectors | update_JxW_values);
    
    // Extractors to real and imaginary parts
    const FEValuesExtractors::Vector E_re(0);
    const FEValuesExtractors::Vector E_im(dim);
    
    std::vector<FEValuesExtractors::Vector> vec(2);
    vec[0] = E_re;
    vec[1] = E_im;
    
    // Local cell storage:
    Vector<double> cell_rhs (dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    
    // Material parameters:
    double current_mur;
    double mur_inv;
    double current_kappa_re;
    double current_kappa_im;
    
    // Neumann storage
    std::vector< Vector<double> > neumann_value_list(n_face_q_points, Vector<double>(fe->n_components()));

    Tensor<1,dim> normal_vector;
    // Storage for both real/imag parts, 0 = re, 1 = im.
    std::vector< Tensor<1,dim> > neumann_value(2);
    std::vector< Tensor<1,dim> > crossproduct_result(2);
    
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_rhs = 0.0;
      current_mur = EquationData::param_mur(cell->material_id());
      mur_inv = 1.0/current_mur;
      current_kappa_re = EquationData::param_kappa_re(cell->material_id());
      current_kappa_im = EquationData::param_kappa_im(cell->material_id());
      // Store coefficients of the real/imaginary blocks:
      std::vector<std::vector<double>> kappa_matrix(2, std::vector<double> (2));
      kappa_matrix[0][0] = current_kappa_re;
      kappa_matrix[0][1] = -current_kappa_im;
      kappa_matrix[1][0] = current_kappa_im;
      kappa_matrix[1][1] = current_kappa_re;
      
      // Loop over faces for neumann condition:
      for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
      {
        fe_face_values.reinit (cell, face_number);
        
        if (cell->face(face_number)->at_boundary()
          &&
          (cell->face(face_number)->boundary_id() == 10))
        {
          // Store values of (mur^-1)*curl E:
          // For this problem, vector value list returns values of H
          // Note that H = i(omega/mu)*curl(E), (mu NOT mur, remember mur = mu/mu0)
          // so (1/mur)*curl(E) = mu_0*H/(i*omega). (1/i = -i)
          // i.e. we use the imaginary part of H for real curl E and real for imag curl E.
          //      and must multiply the imag part by -1.

          boundary_function.curl_value_list(fe_face_values.get_quadrature_points(), neumann_value_list);
          for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
          {
            for (unsigned int d=0; d<dim; ++d)
            {
              neumann_value[0][d] = neumann_value_list[q_point](d);
              neumann_value[1][d] = neumann_value_list[q_point](d+dim);
              normal_vector[d] = fe_face_values.normal_vector(q_point)(d);
            }
            
            for (unsigned int b=0; b<2; ++b)
            {
              cross_product(crossproduct_result[b], normal_vector, neumann_value[b]);
            }
            
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const unsigned int block_index_i = fe->system_to_block_index(i).first;
              
              cell_rhs(i)
              -= mur_inv
              *crossproduct_result[block_index_i]
              *fe_face_values[vec[block_index_i]].value(i,q_point)
              *fe_face_values.JxW(q_point);
            }
          }
        }
      }
      cell_matrix = 0;
      cell->get_dof_indices (local_dof_indices);
      
      // Create the columns of the local matrix which belong to a constrained
      // DoF. These are all that are needed to apply the constraints to the RHS
      fe_values.reinit (cell);
      
      for (unsigned int j=0; j<dofs_per_cell; ++j)
      {
        // Check each column to see if it corresponds to a constrained DoF:
        if ( constraints.is_constrained(local_dof_indices[j]) )
        {
          const unsigned int block_index_j = fe->system_to_block_index(j).first;
          // If yes, cycle through all rows to fill the column:
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int block_index_i = fe->system_to_block_index(i).first;
            
            double curl_part = 0;
            double mass_part = 0;
            if (block_index_i == block_index_j)
            {
              for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              {
                curl_part
                += fe_values[vec[block_index_i]].curl(i, q_point)
                *fe_values[vec[block_index_j]].curl(j, q_point)
                *fe_values.JxW(q_point);
                  
                mass_part
                += fe_values[vec[block_index_i]].value(i, q_point)
                *fe_values[vec[block_index_j]].value(j, q_point)
                *fe_values.JxW(q_point);
              }
              cell_matrix(i,j) += mur_inv*curl_part + kappa_matrix[block_index_i][block_index_j]*mass_part;
            }
            else // off diagonal - curl-curl operator not needed.
            {
              for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              {
                mass_part
                += fe_values[vec[block_index_i]].value(i,q_point)
                *fe_values[vec[block_index_j]].value(j,q_point)
                *fe_values.JxW(q_point);
              }
              cell_matrix(i,j) += kappa_matrix[block_index_i][block_index_j]*mass_part;
            }
          }
        }
      }
      // Use the cell matrix constructed for the constrained DoF columns
      // to add the local contribution to the global RHS:
      constraints.distribute_local_to_global(cell_rhs, local_dof_indices, system_rhs, cell_matrix);      
    }
  }
  
  // Version of assemble_rhs to include the option of a non-zero rhs function
  // This is added as an additional input Function (curlFunctions, etc are automatically valid).
  template <int dim, class DH>
  void EddyCurrent<dim, DH>::assemble_rhs (const DH &dof_handler,
                                           const curlFunction<dim> &boundary_function,
                                           const Function<dim> &rhs_function)
  {
    /*
     * Function to assemble the RHS for a given boundary function.
     * 
     * It first updates the constraints and then computes the RHS.
     */
    
    // Zero the RHS:
    system_rhs = 0;
    
    compute_constraints(dof_handler,
                        boundary_function);
    
    QGauss<dim>  quadrature_formula(quad_order);
    const unsigned int n_q_points = quadrature_formula.size();
    
    QGauss<dim-1> face_quadrature_formula(quad_order);
    const unsigned int n_face_q_points = face_quadrature_formula.size();
    
    const unsigned int dofs_per_cell = fe->dofs_per_cell;

    // Needed to calc the local matrix for distribute_local_to_global
    // Note: only need the columns of the constrained entries.
    FEValues<dim> fe_values (*fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points  |  update_JxW_values);
    
    FEFaceValues<dim> fe_face_values(*fe, face_quadrature_formula,
                                     update_values | update_quadrature_points |
                                     update_normal_vectors | update_JxW_values);
    
    // Extractors to real and imaginary parts
    const FEValuesExtractors::Vector E_re(0);
    const FEValuesExtractors::Vector E_im(dim);
    
    std::vector<FEValuesExtractors::Vector> vec(2);
    vec[0] = E_re;
    vec[1] = E_im;
    
    // Local cell storage:
    Vector<double> cell_rhs (dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    
    // Material parameters:
    double current_mur;
    double mur_inv;
    double current_kappa_re;
    double current_kappa_im;
    
    //RHS storage:
    std::vector<Vector<double> > rhs_value_list(n_q_points, Vector<double>(fe->n_components()));
    std::vector< Tensor<1,dim> > rhs_value(2);

    // Neumann storage
    std::vector< Vector<double> > neumann_value_list(n_face_q_points, Vector<double>(fe->n_components()));

    Tensor<1,dim> normal_vector;
    std::vector< Tensor<1,dim> > neumann_value(2);
    std::vector< Tensor<1,dim> > crossproduct_result(2);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_rhs = 0.0;
      current_mur = EquationData::param_mur(cell->material_id());
      mur_inv = 1.0/current_mur;
      current_kappa_re = EquationData::param_kappa_re(cell->material_id());
      current_kappa_im = EquationData::param_kappa_im(cell->material_id());
      // Store coefficients of the real/imaginary blocks:
      std::vector<std::vector<double>> kappa_matrix(2, std::vector<double> (2));
      kappa_matrix[0][0] = current_kappa_re;
      kappa_matrix[0][1] = -current_kappa_im;
      kappa_matrix[1][0] = current_kappa_im;
      kappa_matrix[1][1] = current_kappa_re;
      
      // Loop over faces for neumann condition:
      for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
      {
        fe_face_values.reinit (cell, face_number);
        
        if (cell->face(face_number)->at_boundary()
          &&
          (cell->face(face_number)->boundary_id() == 10))
        {
          // Store values of (mur^-1)*curl E:
          // For this problem, vector value list returns values of H
          // Note that H = i(omega/mu)*curl(E), (mu NOT mur, remember mur = mu/mu0)
          // so (1/mur)*curl(E) = mu_0*H/(i*omega). (1/i = -i)
          // i.e. we use the imaginary part of H for real curl E and real for imag curl E.
          //      and must multiply the imag part by -1.
          
          // TODO: This is likely to cause problems - May need to separate to a neumann and dirichlet function.
          // For now, we've defined a derived Function class, curlFunction.
          // This seems to get around problems with adding the curl to a Function class.
          
          // NOTE: Have reversed the switching of re and im parts.. was causing issues when using the
          // A formulation properly this needs to be clarified fully and sorted out when I create the definitive version of this
          // class.
          boundary_function.curl_value_list(fe_face_values.get_quadrature_points(), neumann_value_list);
          for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
          {
            for (unsigned int d=0; d<dim; ++d)
            {
              neumann_value[0][d] = neumann_value_list[q_point](d);
              neumann_value[1][d] = neumann_value_list[q_point](d+dim);
              normal_vector[d] = fe_face_values.normal_vector(q_point)(d);
            }
            for (unsigned int b=0; b<2; ++b)
            {
              cross_product(crossproduct_result[b], normal_vector, neumann_value[b]);
            }
            
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const unsigned int block_index_i = fe->system_to_block_index(i).first;
              
              cell_rhs(i)
              -= mur_inv
              *crossproduct_result[block_index_i]
              *fe_face_values[vec[block_index_i]].value(i,q_point)
              *fe_face_values.JxW(q_point);
            }
          }
        }
      }
      
      // Create the columns of the local matrix which belong to a constrained
      // DoF. These are all that are needed to apply the constraints to the RHS
      fe_values.reinit (cell);
      
      // Loop for non-zero right hand side in equation:
      rhs_function.vector_value_list(fe_values.get_quadrature_points(), rhs_value_list);
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          rhs_value[0][d] = rhs_value_list[q_point](d);
          rhs_value[1][d] = rhs_value_list[q_point](d+dim);
        }
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int block_index_i = fe->system_to_block_index(i).first;
          
          cell_rhs(i)
          += rhs_value[block_index_i]
          *fe_values[vec[block_index_i]].value(i,q_point)
          *fe_values.JxW(q_point);
        }
      }
      
      cell_matrix = 0;
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int j=0; j<dofs_per_cell; ++j)
      {
        // Check each column to see if it corresponds to a constrained DoF:
        if ( constraints.is_constrained(local_dof_indices[j]) )
        {
          const unsigned int block_index_j = fe->system_to_block_index(j).first;
          // If yes, cycle through all rows to fill the column:
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int  block_index_i = fe->system_to_block_index(i).first;

            double curl_part = 0;
            double mass_part = 0;
            if (block_index_i == block_index_j)
            {
              for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              {
                curl_part
                += fe_values[vec[block_index_i]].curl(i, q_point)
                *fe_values[vec[block_index_j]].curl(j, q_point)
                *fe_values.JxW(q_point);
                  
                mass_part
                += fe_values[vec[block_index_i]].value(i, q_point)
                *fe_values[vec[block_index_j]].value(j, q_point)
                *fe_values.JxW(q_point);
              }
              cell_matrix(i,j) += mur_inv*curl_part + kappa_matrix[block_index_i][block_index_j]*mass_part;
            }
            else
            {
              for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              {
                mass_part
                += fe_values[vec[block_index_i]].value(i,q_point)
                *fe_values[vec[block_index_j]].value(j,q_point)
                *fe_values.JxW(q_point);
              }
              cell_matrix(i,j) += kappa_matrix[block_index_i][block_index_j]*mass_part;
            }
          }
        }
      }
      // Use the cell matrix constructed for the constrained DoF columns
      // to add the local contribution to the global RHS:
      constraints.distribute_local_to_global(cell_rhs, local_dof_indices, system_rhs, cell_matrix);      
    }
  }
  
  template<int dim, class DH>
  void EddyCurrent<dim, DH>::initialise_solver() // removed:(const bool direct_solver_flag)
  {
    // Initialise the direct solver or preconditioner depending on solver choice.
    // Always use the direct solver for p=0. // TODO update documentation here.
    // Once initialised, the boolean variable initialised will be set to true, so that
    // the initialisation won't be repeated.
    
    
    if (initialised)
    {
      // Nothing to do.
      // TODO: add reset_solver() which deletes the preconditioner and resets the direct solver.
      return;
    }
    else
    {
      if (direct) // use direct solver (LU factorisation of whole matrix).
      {
        direct_solve.initialize(system_matrix);
      }
      // Otherwise:
      // Set up the preconditoner
      // the constructor computes everything required so
      // it can be passed straight to a solve() routine.
      else if (p_order == 0) // use low order preconditioner
      {
        preconditioner = new Preconditioner::EddyCurrentPreconditioner_1x1_lowOrder (system_preconditioner);
      }
      else // p>0 - use low & higher order preconditioner with gradients
      {
        preconditioner = new Preconditioner::EddyCurrentPreconditioner_3x3_lowHighOrderGradients (system_preconditioner);
      }
      initialised = true;
    }
  }
  
  template<int dim, class DH>
  void EddyCurrent<dim, DH>::solve (Vector<double> &output_solution)
  {
    // Solve the linear system for the Eddy Current problem.
    // Can call initialise_solver() before this routine, otherwise
    // it will do it on the first call.
    if (~initialised)
    {
      initialise_solver();
    }
    if (direct)
    {
      // Use if switch enabled.
      // Must have called initialise_solver.   
      direct_solve.vmult (solution, system_rhs);
      constraints.distribute (solution);
    }
    else if (p_order == 0) // use GMRES with low order block preconditioner
    {
      /*GMRES*/        
      SolverControl solver_control (system_matrix.m(),
                                    1e-6*system_rhs.l2_norm(),
                                    true, true); // Add to see residual history
      
      GrowingVectorMemory<BlockVector<double> > vector_memory;
      SolverGMRES<BlockVector<double> >::AdditionalData gmres_data;
      gmres_data.max_n_tmp_vectors = 250;
      gmres_data.right_preconditioning = false;
      gmres_data.force_re_orthogonalization = false;
      
      SolverGMRES<BlockVector<double> > gmres(solver_control, vector_memory,
                                              gmres_data);
      
      gmres.solve(system_matrix,
                  solution,
                  system_rhs,
                  *preconditioner);
      
      // Output iterations to screen
      std::cout << "GMRES Iterations:                "
      << solver_control.last_step()
      << std::endl;
      
      constraints.distribute (solution);
    }
    else // p>0 use GMRES with low & higher order block preconditioner
    {        
      /*GMRES*/        
      SolverControl solver_control (system_matrix.m(),
                                    1e-6*system_rhs.l2_norm(),
                                    true, true); // Add to see residual history
      
      GrowingVectorMemory<BlockVector<double> > vector_memory;
      SolverGMRES<BlockVector<double> >::AdditionalData gmres_data;
      gmres_data.max_n_tmp_vectors = 250;
      gmres_data.right_preconditioning = false;
      gmres_data.force_re_orthogonalization = false;
      
      SolverGMRES<BlockVector<double> > gmres(solver_control, vector_memory,
                                              gmres_data);
      
      gmres.solve(system_matrix,
                  solution,
                  system_rhs,
                  *preconditioner);
      
      // Output iterations to screen
      std::cout << "GMRES Iterations:                "
      << solver_control.last_step()
      << std::endl;
      
      constraints.distribute (solution);
    }
  
    // TODO: may not be able to pass to the output solution like this:
    // seems to be ok. for now.
    output_solution.reinit(solution.size());
    output_solution = solution;
  }
  // Template instantiation
  template class EddyCurrent<3, DoFHandler<3>>;
  // END CLASS EDDYCURRENT
}
// END namespace ForwardSolver
