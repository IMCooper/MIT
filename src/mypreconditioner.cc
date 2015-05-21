#include <mypreconditioner.h>

using namespace dealii;

namespace Preconditioner
{
  // base class:
  EddyCurrentPreconditionerBase::EddyCurrentPreconditionerBase (const BlockSparseMatrix<double> &A)
  :
  precon_matrix (&A)
  {
  }
  void EddyCurrentPreconditionerBase::vmult (BlockVector<double>       &dst,
                                             const BlockVector<double> &src) const
  {
    precon_matrix->vmult(dst,src);
  }
  void EddyCurrentPreconditionerBase::Tvmult (BlockVector<double>       &dst,
                                             const BlockVector<double> &src) const
  {
    precon_matrix->Tvmult(dst,src);
  }
  // END base class
  // EddyCurrentPreconditioner_1x1_lowOrder
  EddyCurrentPreconditioner_1x1_lowOrder::EddyCurrentPreconditioner_1x1_lowOrder (const BlockSparseMatrix<double> &A)
  :
  EddyCurrentPreconditionerBase(A)
//   precon_matrix (&A)
//   tmp0 (A.block(0,0).m())
  {
    tmp0.reinit(precon_matrix->block(0,0).m());
    // Block 0:
    //direct:
    block0.initialize(precon_matrix->block(0,0));

    // ILU:
//     block0_data.extra_off_diagonals = PreconditionerData::extra_off_diagonals;
//     block0_data.strengthen_diagonal = PreconditionerData::strengthen_diagonal;
//     block0.initialize(precon_matrix->block(0,0),
//                       block0_data);
  }
  void EddyCurrentPreconditioner_1x1_lowOrder::vmult (BlockVector<double>       &dst,
                                                      const BlockVector<double> &src) const
  {
    block0.vmult (tmp0, src.block(0));
    dst.block(0)=tmp0;
  }
  void EddyCurrentPreconditioner_1x1_lowOrder::Tvmult (BlockVector<double>       &dst,
                                                       const BlockVector<double> &src) const
  {
    block0.Tvmult (tmp0, src.block(0));
    dst.block(0)=tmp0;
  }
  // END EddyCurrentPreconditioner_1x1_lowOrder

  // EddyCurrentPreconditioner_2x2_lowAndHighOrder
  EddyCurrentPreconditioner_2x2_lowHighOrder::EddyCurrentPreconditioner_2x2_lowHighOrder (const BlockSparseMatrix<double> &A)
  :
  EddyCurrentPreconditionerBase(A)
  {
    tmp0.reinit(precon_matrix->block(0,0).m());
    tmp1.reinit(precon_matrix->block(1,1).m());
    // Block 0:
    //direct:
    block0.initialize(precon_matrix->block(0,0));

    // ILU:
//     block0_data.extra_off_diagonals = PreconditionerData::extra_off_diagonals;
//     block0_data.strengthen_diagonal = PreconditionerData::strengthen_diagonal;
//     block0.initialize(precon_matrix->block(0,0),
//                       block0_data);

    // SOR:
//     block0.initialize(precon_matrix->block(0,0),1e6);

    //Block 1:
    // Options for ILU:
    block1_data.extra_off_diagonals = PreconditionerData::extra_off_diagonals;
    block1_data.strengthen_diagonal = PreconditionerData::strengthen_diagonal;

    block1.initialize(precon_matrix->block(1,1),
                      block1_data);
  }
  void EddyCurrentPreconditioner_2x2_lowHighOrder::vmult (BlockVector<double>       &dst,
                                         const BlockVector<double> &src) const
  {
    block0.vmult (tmp0, src.block(0));
    block1.vmult (tmp1, src.block(1));

    dst.block(0)=tmp0;
    dst.block(1)=tmp1;
  }
  void EddyCurrentPreconditioner_2x2_lowHighOrder::Tvmult (BlockVector<double>       &dst,
                                          const BlockVector<double> &src) const
  {
    block0.Tvmult (tmp0, src.block(0));
    block1.Tvmult (tmp1, src.block(1));

    dst.block(0)=tmp0;
    dst.block(1)=tmp1;

  }
  // END EddyCurrentPreconditioner_2x2_lowHighOrder
  
  // EddyCurrentPreconditioner_3x3_lowHighOrderGradients
    EddyCurrentPreconditioner_3x3_lowHighOrderGradients::EddyCurrentPreconditioner_3x3_lowHighOrderGradients (const BlockSparseMatrix<double> &A)
  :
  EddyCurrentPreconditionerBase(A)
  {
    tmp0.reinit(precon_matrix->block(0,0).m());
    tmp1.reinit(precon_matrix->block(1,1).m());
    tmp2.reinit(precon_matrix->block(2,2).m());
    // Block 0:
    //direct:
    block0.initialize(precon_matrix->block(0,0));
    
    // ILU:
//     block0_data.extra_off_diagonals = PreconditionerData::extra_off_diagonals;
//     block0_data.strengthen_diagonal = PreconditionerData::strengthen_diagonal;
//     block0.initialize(precon_matrix->block(0,0),
//                       block0_data);

    //Block 1:
    // Options for ILU:
//     block1_data.extra_off_diagonals = PreconditionerData::extra_off_diagonals;
//     block1_data.strengthen_diagonal = PreconditionerData::strengthen_diagonal;

//     block1.initialize(precon_matrix->block(1,1),
//                       block1_data);
    block1.initialize(precon_matrix->block(1,1));
    
    //Block 2:
    // Options for ILU:
//     block2_data.extra_off_diagonals = PreconditionerData::extra_off_diagonals;
//     block2_data.strengthen_diagonal = PreconditionerData::strengthen_diagonal;

//     block2.initialize(precon_matrix->block(2,2),
//                       block2_data);
    block2.initialize(precon_matrix->block(2,2));

  }
  void EddyCurrentPreconditioner_3x3_lowHighOrderGradients::vmult (BlockVector<double>       &dst,
                                         const BlockVector<double> &src) const
  {
    block0.vmult (tmp0, src.block(0));
    block1.vmult (tmp1, src.block(1));
    block2.vmult (tmp2, src.block(2));

    dst.block(0)=tmp0;
    dst.block(1)=tmp1;
    dst.block(2)=tmp2;
  }
  void EddyCurrentPreconditioner_3x3_lowHighOrderGradients::Tvmult (BlockVector<double>       &dst,
                                          const BlockVector<double> &src) const
  {
    block0.Tvmult (tmp0, src.block(0));
    block1.Tvmult (tmp1, src.block(1));
    block2.Tvmult (tmp2, src.block(2));

    dst.block(0)=tmp0;
    dst.block(1)=tmp1;
    dst.block(2)=tmp2;
  }
  // END EddyCurrentPreconditioner_3x3_lowHighOrderGradients
}
// END namespace Preconditioner