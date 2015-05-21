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
  // base class for all eddy current preconditioners:
  class EddyCurrentPreconditionerBase : public Subscriptor
  {
  public:
    EddyCurrentPreconditionerBase (const BlockSparseMatrix<double> &A);
    
    virtual void vmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
    virtual void Tvmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
  protected:
    const SmartPointer<const BlockSparseMatrix<double> > precon_matrix;
  };

  // Preconditioner for iterative solver, made up of a single block.
  // This essentially assumes the use of p=0.
  class EddyCurrentPreconditioner_1x1_lowOrder : public virtual EddyCurrentPreconditionerBase
  {
  public:
    EddyCurrentPreconditioner_1x1_lowOrder (const BlockSparseMatrix<double> &A);

    virtual void vmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
    virtual void Tvmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
  private:
    mutable Vector<double> tmp0;

    // Various options for the single block.
//     PreconditionIdentity block0;
//     PreconditionJacobi<> block0;
//     SparseILU<double> block0;
//     SparseILU<double>::AdditionalData block0_data;
//     PreconditionSOR<SparseMatrix<double> > block0;
    SparseDirectUMFPACK block0;
  };

  // Preconditioning for use with an FE_Nedelec-based 2x2 block matrix
  // made up of:
  // block 0: lowest
  // block 1: higher order edges/faces/cells.
  class EddyCurrentPreconditioner_2x2_lowHighOrder : public EddyCurrentPreconditionerBase
  {
  public:
    EddyCurrentPreconditioner_2x2_lowHighOrder (const BlockSparseMatrix<double> &A);

    void vmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
    void Tvmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
  private:
//     const SmartPointer<const BlockSparseMatrix<double> > precon_matrix;
    mutable Vector<double> tmp0, tmp1;

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
  };

  // Preconditioning for use with an FE_Nedelec-based 3x3 block matrix
  // made up of:
  // block 0: lowest
  // block 1: gradient-based higher order edges/faces/cells
  // block 2: non-gradient-based higher order edges/faces/cells
  // TODO: Make use of a nested CG solver for the higher order blocks.
  class EddyCurrentPreconditioner_3x3_lowHighOrderGradients : public EddyCurrentPreconditionerBase
  {
  public:
    EddyCurrentPreconditioner_3x3_lowHighOrderGradients (const BlockSparseMatrix<double> &A);

    void vmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
    void Tvmult (BlockVector<double> &dst, const BlockVector<double> &src) const;
  private:
//     const SmartPointer<const BlockSparseMatrix<double> > precon_matrix;
    mutable Vector<double> tmp0, tmp1, tmp2;

    // Block 0: lowest order:
//     PreconditionIdentity block0;
//     PreconditionJacobi<> block0;
//     SparseILU<double> block0;
//     SparseILU<double>::AdditionalData block0_data;
//     PreconditionSOR<SparseMatrix<double> > block0;
    SparseDirectUMFPACK block0;

    // Block 1: higher order gradients:
//     PreconditionIdentity block1;
//     PreconditionJacobi<> block1;
    SparseDirectUMFPACK block1;
//     PreconditionSOR<SparseMatrix<double> > block1;
//     SparseILU<double> block1;
//     SparseILU<double>::AdditionalData block1_data;

    // Block 2: higher order non-gradients:
//     PreconditionIdentity block2;
//     PreconditionJacobi<> block2;
    SparseDirectUMFPACK block2;
//     PreconditionSOR<SparseMatrix<double> > block2;
//     SparseILU<double> block2;
//     SparseILU<double>::AdditionalData block2_data;
  };
}
#endif