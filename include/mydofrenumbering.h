// deal.II includes:
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>


using namespace dealii;

#ifndef MYDOFRENUMBERING_H
#define MYDOFRENUMBERING_H
namespace MyDoFRenumbering
{
  template <int dim, class DH>
  void by_dimension (DH &dof_handler,
                     std::vector<unsigned int> &dof_counts);
}
#endif