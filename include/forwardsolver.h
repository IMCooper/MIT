// deal.II includes:
#include <deal.II/base/quadrature.h>
#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/sparse_decomposition.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

// std includes:

// My includes:
#include <all_data.h>
//#include <backgroundfield.h>
#include <curlfunction.h>
#include <mydofrenumbering.h>
#include <mypreconditioner.h>

#include <omp.h>
#include <Config.h>
#include <Timer.h>
using namespace dealii;

#ifndef FORWARDSOLVER_H
#define FORWARDSOLVER_H
namespace ForwardSolver
{
  using namespace dealii;
  // TODO:
  // At the moment this is a general vector-wave equation solver. ?? It may be a good idea to rename it ???
  //
  // It can be made to be an eddy current solver with the correct choice of material parameters.
  //
  // Care must be taken with the choice of material parameters, boundary conditions and RHS (currently enforced to be zero)
  // in that they must all be chosen to be consistent with each other.
  //
  // i.e. if we use an A-based formulation:
  // curl(1/mu_r*(curl(A)) + kappa*mu_0*A = mu_0*J_s
  //
  // then the A and curl(A) for setting the BCs must be chosen to reflect this
  //
  // similarly, if solving for E:
  // curl(1/mu_r*(curl(E)) + kappa*mu_0*E = i*omega*mu_0*J_s
  //
  // then E and curl(E) must be set appropriately.
  //
  //
  // Also, in the case of neumann BCs, they must be applied using the curlFunction<dim>
  // and then the contribution is n x mu_{r} x curl(F) = ... 
  // where curl(F) is defined through the derived curlFunction<dim> passed to the solver.
  // note F is either A or E here.
  //
  // Finally, care in the choice of regularisation must be taken. Where kappa=0 then we should choose a small regularisation
  // to stop the system becoming singular, so that kappa is set to be kappa_re = reg*mu0, kappa_im = 0.
  // 
  
  template <int dim, class DH>
  class EddyCurrent
  {
    // This class is for the solution of the eddy current approximation of the time-harmonic maxwell equations.
    //  
    // Input:
    // - a dof_handler which is attached to an finite element and triangulation with a set of material IDs and BC IDs.
    // - A list of material parameter values corresponding to these IDs.
    //  
    // Functions:
    // - compute constraints on the linear system according to the BCs and an input function for the BCs.
    // - compute the stiffness matrix associated with the PDE
    // - compute the RHS vector associated with the PDE
    // - solve the linear system.
    // 
    // Output:
    // - a solution vector for the coeffs of the FE solution of the PDE
  public:
    EddyCurrent (DH &dof_handler,
                 const FiniteElement<dim> &fe,
                 imc::Config* aConfig, 
                 const bool direct_flag = false);
    EddyCurrent (const Mapping<dim,dim> &mapping_in,
                 DH &dof_handler,
                 const FiniteElement<dim> &fe,
                 imc::Config* aConfig, 
                 const bool direct_flag = false);
    
    ~EddyCurrent ();
    
    void assemble_matrices (const DH &dof_handler);
    void assemble_matrices_OMP (const DH &dof_handler);
    // zero rhs version:
    void assemble_rhs (const DH &dof_handler,
                       const curlFunction<dim> &boundary_function);
    void assemble_rhs_OMP (const DH &dof_handler,
                       const curlFunction<dim> &boundary_function);

    
    // non-zero rhs version:
    void assemble_rhs (const DH &dof_handler,
                       const curlFunction<dim> &boundary_function,
                       const Function<dim> &rhs_function);
    // non-zero rhs version:
    void assemble_rhs_OMP (const DH &dof_handler,
                       const curlFunction<dim> &boundary_function,
                       const Function<dim> &rhs_function);
    
    void initialise_solver();
    void solve (Vector<double> &output_solution,
                unsigned int &n_iterations);
    
    
  private:
    void constructor_setup(DH &dof_handler,
                           const bool direct_flag);
    
    void compute_constraints (const DH &dof_handler,
                              const curlFunction<dim> &boundary_function);
    
    // FE storage:
    const SmartPointer< const FiniteElement<dim> > fe; // TODO: check if ok, or should be protected??
    unsigned int p_order;
    unsigned int quad_order;
    
    imc::Config* config;    

    // Mapping storage:
    SmartPointer< const Mapping<dim> > mapping;
    
    // Block sparse storage:
    BlockSparsityPattern sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
    BlockVector<double> solution;
    BlockVector<double> system_rhs;
    
    // Dof ordering storage:
    unsigned int n_lowest_order_dofs;
    unsigned int n_higher_order_dofs;
    
    unsigned int n_higher_order_edge_gradients_dofs;
    unsigned int n_higher_order_face_gradients_dofs;
    unsigned int n_higher_order_cell_gradients_dofs;
    unsigned int n_higher_order_face_nongradients_dofs;
    unsigned int n_higher_order_cell_nongradients_dofs;
    
    unsigned int n_higher_order_gradient_dofs;
    unsigned int n_higher_order_non_gradient_dofs;
    
    // Preconditioner:
    BlockSparseMatrix<double> system_preconditioner;
    // TODO: would this be better as SmartPointer<Preconditioner::EddyCurrentPreconditionerBase> preconditioner; ??
    Preconditioner::EddyCurrentPreconditionerBase* preconditioner;
    
    // Constraints:
    ConstraintMatrix constraints;
    
    // Linear Algebra Storage:
    bool initialised = false;
    bool direct = false;
    SparseDirectUMFPACK direct_solve;
    
    unsigned int omp_threads;
    unsigned int omp_chunkSize;
    unsigned int omp_collapseLoops;
  
    // protected: // removed, used to contain smart pointer for fe.
  };  
}
#endif
