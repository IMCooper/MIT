#include <mypreconditioner.h>

using namespace dealii;

namespace Preconditioner
{
  // HIGHER ORDER VERSION:
  EddyCurrentPreconditioner::EddyCurrentPreconditioner (const BlockSparseMatrix<double> &A)
  :
  precon_matrix (&A),
  tmp0 (A.block(0,0).m()),
  tmp1 (A.block(1,1).m())
//   tmp2 (A.block(2,2).m()),
//   tmp3 (A.block(3,3).m())
  {
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
    
    // SOR:
//     block1.initialize(precon_matrix->block(1,1),1e6);

//     block2.initialize(precon_matrix->block(2,2));

//     block3.initialize(precon_matrix->block(3,3));

  }
  void EddyCurrentPreconditioner::vmult (BlockVector<double>       &dst,
                                         const BlockVector<double> &src) const
  {
    block0.vmult (tmp0, src.block(0));
    block1.vmult (tmp1, src.block(1));
//     block2.vmult (tmp2, src.block(2));
//     block3.vmult (tmp3, src.block(3));
    dst.block(0)=tmp0;
    dst.block(1)=tmp1;
//     dst.block(2)=tmp2;
//     dst.block(3)=tmp3;
    
  }
  void EddyCurrentPreconditioner::Tvmult (BlockVector<double>       &dst,
                                          const BlockVector<double> &src) const
  {
    block0.Tvmult (tmp0, src.block(0));
    block1.Tvmult (tmp1, src.block(1));
//     block2.Tvmult (tmp2, src.block(2));
//     block3.Tvmult (tmp3, src.block(3));
    dst.block(0)=tmp0;
    dst.block(1)=tmp1;
//     dst.block(2)=tmp2;
//     dst.block(3)=tmp3;
  }
  // END EDDYCURRENTPRECONDITIONER
  
  
  // EDDYCURRENTPRECONDITIONER_LOW_ORDER
  EddyCurrentPreconditioner_low_order::EddyCurrentPreconditioner_low_order (const BlockSparseMatrix<double> &A)
  :
  precon_matrix (&A),
  tmp0 (A.block(0,0).m())
  {
    // Block 0:
    //direct:
    block0.initialize(precon_matrix->block(0,0));
    
    // ILU:
//     block0_data.extra_off_diagonals = PreconditionerData::extra_off_diagonals;
//     block0_data.strengthen_diagonal = PreconditionerData::strengthen_diagonal;
//     block0.initialize(precon_matrix->block(0,0),
//                       block0_data);
    
    // SOR:
//     block0.initialize(precon_matrix->block(0,0),(1.0/6.28e7));
    


  }
  void EddyCurrentPreconditioner_low_order::vmult (BlockVector<double>       &dst,
                                         const BlockVector<double> &src) const
  {
    block0.vmult (tmp0, src.block(0));
    dst.block(0)=tmp0;    
  }
  void EddyCurrentPreconditioner_low_order::Tvmult (BlockVector<double>       &dst,
                                          const BlockVector<double> &src) const
  {
    block0.Tvmult (tmp0, src.block(0));
    dst.block(0)=tmp0;
  }
  // END EDDYCURRENTPRECONDITIONER_LOW_ORDER
}
// END namespace Preconditioner