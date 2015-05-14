#include <mydofrenumbering.h>



namespace MyDoFRenumbering
{
  using namespace dealii;
  namespace internals
  {
    // internal function to do the work:
    template <int dim, class DH>
    types::global_dof_index
    compute_by_dimension (std::vector<types::global_dof_index> &renumbering,
                          const DH &dof_handler,
                          std::vector<unsigned int> &dof_counts)
    {
      // Splits DoFs for an FENedelec<dim> into lower order edges and higher order edges/faces/cells
      // 
      // For each copy of each base element, we know there are 12 lowest order
      // edge DoFs on the reference element
      //  
      // Loop through each base element and each multiplicity and grab the system index of each.
      // Can then compare this index to where we know the lowest order DoF is and reorder. 
      
      const FiniteElement<dim> &fe = dof_handler.get_fe ();
      const unsigned int superdegree = fe.degree;
      const unsigned int degree = superdegree - 1;
      const unsigned int n_global_dofs = dof_handler.n_dofs();
      const unsigned int lines_per_cell = GeometryInfo<dim>::lines_per_cell;
      const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
      
      const unsigned int dofs_per_line = fe.dofs_per_line;
      const unsigned int dofs_per_quad = fe.dofs_per_quad;
      const unsigned int dofs_per_hex = fe.dofs_per_hex;
      
      const unsigned int dofs_per_cell = fe.dofs_per_cell;
      
      const unsigned int line_dofs_per_cell = lines_per_cell*dofs_per_line;
      const unsigned int face_dofs_per_cell = faces_per_cell*dofs_per_quad;
      
      
      std::vector<bool> global_dofs_seen (n_global_dofs);
      
      /*
       * Store a temporary array to store all renumberings
       * as well as an array containing a count of how many DoFs have
       * been reordered within each
       * we'll then combine these lists to create the final cell_reordering
       * with a loop over edges, faces and then cells, with the lowest order
       * being first within each.
       */
      
      
      
      for (unsigned int i = 0; i < n_global_dofs; ++i)
      {
        global_dofs_seen[i] = false;
      }
      
      std::pair<unsigned int, unsigned int> base_indices (0,0);
      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
      /*
       * Idea is:
       * We know that for an FESystem of Nedelec elements we have a local ordering
       * based on 2 (in 2D) or 3 (in 3D) loops:
       * loop 1:
       * base elements
       * -> element copy
       *    -> edges
       *       -> basis polynomial orders
       * 
       * loop 2:
       * base elements
       * -> element copy
       *    -> faces
       *       -> basis polynomial orders
       * 
       * loop 3:
       * base elements
       * -> element copy
       *    -> cells
       *       -> basis polynomial orders
       * 
       * So if we want to group the DoFs by base element, element copy and polynomial order on edges/faces/cells
       * then we can exploit this.
       * 
       * We also know that on a local cell, the DoFs are grouped by:
       * edges then
       * faces then
       * cells
       * 
       * For now we only want to find the lowest order edge DoFs.
       * 
       * within a given copy of a given base element (i.e. we can say base_indices = (base_element, element_copy))
       * we can look at the local cell index, restricting ourselves to edge DoFs since we know they come before faces/cells,
       * (via base_cell_index = fe.system_to_base_index(local_dof_index).second)
       * and know that a lowest order edge function satifies mod(base_cell_index,fe.degree) = 0
       * 
       * similarly, the 1st order (assuming fe.degree>1) edge function satifies
       * mod(base_cell_index,fe.degree) = 1, and so on.
       * 
       * Note: It is not quite so simple for face and cell DoFs, since the number of them as the degree increases
       *       does not scale linearly.
       * 
       */
      
      /*
       * Temporary arrays to store each type of DoF 
       * (lowest order and higher order edge/face/cell).
       *
       * Will combine these into the final reordering afterwards.
       */
      /* eventually want to have reorderings by both edge/face/cell and by polynomial degree within each:
       * 
       * std::vector<std::vector<types::global_dof_index>> edge_reorderings (superdegree, std::vector<types::global_dof_index> (n_global_dofs));
       * std::vector<types::global_dof_index> edge_counts (superdegree);
       * std::vector<std::vector<types::global_dof_index>> face_reorderings (superdegree, std::vector<types::global_dof_index> (n_global_dofs));
       * std::vector<types::global_dof_index> face_counts (superdegree-1);
       * std::vector<std::vector<types::global_dof_index>> cell_reorderings (superdegree, std::vector<types::global_dof_index> (n_global_dofs));
       * std::vector<types::global_dof_index> cell_counts (superdegree-1);
       */
      std::vector<types::global_dof_index> lowest_order_renumbering(n_global_dofs);
      std::vector<types::global_dof_index> higher_order_edge_renumbering(n_global_dofs);
      std::vector<types::global_dof_index> higher_order_face_renumbering(n_global_dofs);
      std::vector<types::global_dof_index> higher_order_cell_renumbering(n_global_dofs);
      
      
      std::vector<bool> lowest_order_dof(n_global_dofs);
      std::vector<bool> higher_order_edge_dof(n_global_dofs);
      std::vector<bool> higher_order_face_dof(n_global_dofs);
      std::vector<bool> higher_order_cell_dof(n_global_dofs);
      
      for (unsigned int dof = 0; dof < n_global_dofs; ++dof)
      {
        lowest_order_dof[dof] = false;
        higher_order_edge_dof[dof] = false;
        higher_order_face_dof[dof] = false;
        higher_order_cell_dof[dof] = false;
      }
      
      types::global_dof_index lowest_order_dof_count = 0;
      types::global_dof_index higher_order_edge_dof_count = 0;
      types::global_dof_index higher_order_face_dof_count = 0;
      types::global_dof_index higher_order_cell_dof_count = 0;
      
      for (unsigned int base_element = 0; base_element < fe.n_base_elements (); ++base_element)
      {
        base_indices.first = base_element;
        for (unsigned int element_copy = 0; element_copy < fe.element_multiplicity (base_element); ++element_copy)
        {
          base_indices.second = element_copy;
          typename DoFHandler<dim>::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();
          for (; cell!=endc; ++cell)
          {
            cell->get_dof_indices (local_dof_indices);
            std::vector<bool> local_dofs_seen (dofs_per_cell);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              local_dofs_seen[i] = false;
            }
            /* 
             * Loop over all edge DoFs.
             * Note that edge dofs are all grouped together
             * so we can separate the lowest orders from
             * the higher orders and not worry about face/cell DoFs.
             */
            for (unsigned int local_dof = 0; local_dof < line_dofs_per_cell; ++local_dof)
            {
              // Lowest order:
              if ( (fe.system_to_base_index(local_dof).first == base_indices)
                && (fe.system_to_base_index(local_dof).second % superdegree == 0)
                && ( !local_dofs_seen[local_dof] )
                && ( !global_dofs_seen[local_dof_indices[local_dof]] ) )
              {
                local_dofs_seen[local_dof] = true;
                global_dofs_seen[local_dof_indices[local_dof]] = true;
                lowest_order_renumbering[local_dof_indices[local_dof]] = lowest_order_dof_count;
                lowest_order_dof[local_dof_indices[local_dof]] = true;
                ++lowest_order_dof_count;               
              }
              else if ( (fe.system_to_base_index(local_dof).first == base_indices)
                && (fe.system_to_base_index(local_dof).second % superdegree > 0)
                && ( !local_dofs_seen[local_dof] )
                && ( !global_dofs_seen[local_dof_indices[local_dof]] ) )
              {
                local_dofs_seen[local_dof] = true;
                global_dofs_seen[local_dof_indices[local_dof]] = true;
                higher_order_edge_renumbering[local_dof_indices[local_dof]] = higher_order_edge_dof_count;
                higher_order_edge_dof[local_dof_indices[local_dof]] = true;
                ++higher_order_edge_dof_count;
              }
            }
            /* 
             * Loop over all face DoFs, skipping the edge DoFs, stopping
             * before we reach the interior cell DoFs.
             * Since all face DoFs are already higher order, and for
             * now we're not reordering by polynomial degree, there is no
             * need to do much else.
             * We simply group these DoFs using their common base element and copy.
             */
            for (unsigned int local_dof = line_dofs_per_cell;
                 local_dof < line_dofs_per_cell + face_dofs_per_cell; ++local_dof)
            {
              if ( (fe.system_to_base_index(local_dof).first == base_indices)
                && ( !local_dofs_seen[local_dof] )
                && ( !global_dofs_seen[local_dof_indices[local_dof]] ) )
              {
                local_dofs_seen[local_dof] = true;
                global_dofs_seen[local_dof_indices[local_dof]] = true;
                higher_order_face_renumbering[local_dof_indices[local_dof]] = higher_order_face_dof_count;
                higher_order_face_dof[local_dof_indices[local_dof]] = true;
                ++higher_order_face_dof_count;
              }
            }
            /*
             * Loop over all cell DoFs, skipping the edge and face DoFs,
             * stopping at the last local DoF
             * Since all cell DoFs are already higher order, and for
             * now we're not reordering by polynomial degree, there is no
             * need to do much else.
             * We simply group these DoFs using their common base element and copy.
             */
            for (unsigned int local_dof = line_dofs_per_cell + face_dofs_per_cell;
                 local_dof < dofs_per_cell; ++local_dof)
            {
              if ( (fe.system_to_base_index(local_dof).first == base_indices)
                && ( !local_dofs_seen[local_dof] )
                && ( !global_dofs_seen[local_dof_indices[local_dof]] ) )
              {
                local_dofs_seen[local_dof] = true;
                global_dofs_seen[local_dof_indices[local_dof]] = true;
                higher_order_cell_renumbering[local_dof_indices[local_dof]] = higher_order_cell_dof_count;
                higher_order_cell_dof[local_dof_indices[local_dof]] = true;
                ++higher_order_cell_dof_count;
              }
            }
          }
        }
      }
      /* Fill the array with counts of each type of DoF */
      dof_counts[0] = lowest_order_dof_count;
      dof_counts[1] = higher_order_edge_dof_count;
      dof_counts[2] = higher_order_face_dof_count;
      dof_counts[3] = higher_order_cell_dof_count;
      /* 
       * Now have reordering for lowest and higher order (split into edges/face/cells) dofs.
       * Now want to combine them into a single reordering in this order:
       *  - lowest order first
       *  - higher order edges
       *  - higher order faces
       *  - higher order cells
       */
      //    renumbering = lowest_order_renumbering;
      // Higher order edges:
      for (unsigned int dof = 0; dof < n_global_dofs; ++dof)
      {
        if (lowest_order_dof[dof])
        {
          renumbering[dof] = lowest_order_renumbering[dof];
        }
        else if (higher_order_edge_dof[dof])
        {
          renumbering[dof] = higher_order_edge_renumbering[dof]
          + lowest_order_dof_count;
        }
        else if (higher_order_face_dof[dof])
        {
          renumbering[dof] = higher_order_face_renumbering[dof]
          + lowest_order_dof_count
          + higher_order_edge_dof_count;
        }
        else if (higher_order_cell_dof[dof])
        {
          renumbering[dof] = higher_order_cell_renumbering[dof]
          + lowest_order_dof_count
          + higher_order_edge_dof_count
          + higher_order_face_dof_count;
        }
      }
      /* Check if we've now reordered all DoFs */
      unsigned int total_counted_dofs;
      total_counted_dofs = lowest_order_dof_count
      + higher_order_edge_dof_count
      + higher_order_face_dof_count
      + higher_order_cell_dof_count;
      if (total_counted_dofs == n_global_dofs)
      {
        return total_counted_dofs;
      }
      // If not, then add remaining DoFs in the order they come:
      //     std::cout << " DoFRenumbering::by_dimension: some DoFs still remain..." << std::endl;
      for (unsigned int dof = 0; dof < n_global_dofs; ++dof)
      {
        if (!global_dofs_seen[dof])
        {
          renumbering[dof] = total_counted_dofs;
          ++total_counted_dofs;
          global_dofs_seen[dof] = true;
        }
      }
      return total_counted_dofs;
    }
  }
  // Short function to call the internal function:
  template <int dim, class DH>
  void by_dimension (DH &dof_handler,
                     std::vector<unsigned int> &dof_counts)
  {
    std::vector<types::global_dof_index> renumbering(dof_handler.n_dofs());
    
    types::global_dof_index result = internals::
                                     compute_by_dimension<dim, DH>(renumbering,
                                                                   dof_handler,
                                                                   dof_counts);
    if (result == 0)
    {
      return;
    }
    dof_handler.renumber_dofs (renumbering);    
  }
  
  template
  void by_dimension<3>(DoFHandler<3> &,
                      std::vector<unsigned int> &);  
}