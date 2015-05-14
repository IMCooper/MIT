// deal.II includes:
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/sparse_decomposition.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

// std includes:

// My includes:
#include <all_data.h>

using namespace dealii;

#ifndef MYPRECONDITIONER_H
#define MYPRECONDITIONER_H
namespace Preconditioner
{
  // Preconditioning for use with an FE_Nedelec-based block matrix made up of lowest/higher order edges/faces/cells
  // If we further split the blocks by real/imaginary parts then we'd need 8 blocks in total.
  using namespace dealii;
  class EddyCurrentPreconditioner : public Subscriptor
  {
  public:
    EddyCurrentPreconditioner (const BlockSparseMatrix<double> &A);
    
    void vmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
    void Tvmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
  private:
    const SmartPointer<const BlockSparseMatrix<double> > precon_matrix;
    mutable Vector<double> tmp0, tmp1;//, tmp2, tmp3;
    
    // Block 0: lowest order edges.    
//     PreconditionIdentity block0;
//     PreconditionJacobi<> block0;
//     SparseILU<double> block0;
//     SparseILU<double>::AdditionalData block0_data;
//     PreconditionSOR<SparseMatrix<double> > block0;
    SparseDirectUMFPACK block0;    
    
    // Block 1: higher order edges/faces/cells (or just edges if Block2/3 enabled:
//     PreconditionIdentity block1;
//      PreconditionJacobi<> block1;
    SparseILU<double> block1;
    SparseILU<double>::AdditionalData block1_data;
//     PreconditionSOR<SparseMatrix<double> > block1;
//    SparseDirectUMFPACK block1;
  
    // Block 2: higher order faces:
//     PreconditionIdentity block2;
//      PreconditionJacobi<> block2;
//     SparseILU<double> block2;
//    SparseDirectUMFPACK block2;
    
    // Block 3: higher order cells:
//     PreconditionIdentity block3;
//      PreconditionJacobi<> block3;
//     SparseILU<double> block3;
//    SparseDirectUMFPACK block3;
  };
  
    
  // Preconditioner for iterative solver with 1 block (i.e. p=0 case):
  class EddyCurrentPreconditioner_low_order : public Subscriptor
  {
  public:
    EddyCurrentPreconditioner_low_order (const BlockSparseMatrix<double> &A);
    
    void vmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
    void Tvmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
  private:
    const SmartPointer<const BlockSparseMatrix<double> > precon_matrix;
    mutable Vector<double> tmp0;
    
//     PreconditionIdentity block0;
//     PreconditionJacobi<> block0;
//     SparseILU<double> block0;
//     SparseILU<double>::AdditionalData block0_data;
//     PreconditionSOR<SparseMatrix<double> > block0;
    SparseDirectUMFPACK block0;
  };
  
} // END namespace Preconditioner
#endif