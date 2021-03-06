#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <all_data.h>
#include <backgroundfield.h>
#include <curlfunction.h>
#include <forwardsolver.h>
#include <inputtools.h>
#include <mydofrenumbering.h>
#include <mypreconditioner.h>
#include <myvectortools.h>
#include <outputtools.h>

using namespace dealii;

namespace cubeBenchmark
{
  
  template <int dim>
  class cubeBenchmark
  {
  public:
    cubeBenchmark (const unsigned int order);
    ~cubeBenchmark ();
    void run(std::string input_filename,
             std::string output_filename);
  private:
    Triangulation<dim> tria;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;

    unsigned int p_order;
    
    void initialise_materials();
    void process_mesh(bool neuman_flag);
  };
  
  template <int dim>
  cubeBenchmark<dim>::cubeBenchmark(const unsigned int order)
  :
  fe (FE_Nedelec<dim>(order), 2),
  dof_handler (tria),
  p_order(order)
  {
  }
  
  template <int dim>
  cubeBenchmark<dim>::~cubeBenchmark ()
  {
    dof_handler.clear ();  
  }
  
  template <int dim>
  void cubeBenchmark<dim>::initialise_materials()
  {
    EquationData::param_mur.reinit(2);
    EquationData::param_mur(0) = EquationData::param_mur_background;       
    EquationData::param_mur(1) = EquationData::param_mur_conducting;
    
    EquationData::param_sigma.reinit(2);
    EquationData::param_sigma(0) = EquationData::param_sigma_background;
    EquationData::param_sigma(1) = EquationData::param_sigma_conducting;
    
    EquationData::param_epsilon.reinit(2);
    EquationData::param_epsilon(0) = EquationData::param_epsilon_background;
    EquationData::param_epsilon(1) = EquationData::param_epsilon_conducting;
    
    
    // kappa = -omega^2*epr + i*omega*sigma;
    // i.e. kappa_re = -omega^2
    //      kappa_im = omega*sigma
    EquationData::param_kappa_re.reinit(EquationData::param_mur.size());
    EquationData::param_kappa_im.reinit(EquationData::param_mur.size());
    for (unsigned int i=0;i<EquationData::param_mur.size();i++) // note mur and kappa must have same size:
    {
      if (EquationData::param_sigma(i) > 0)
      {
        EquationData::param_kappa_re(i) = 0.0;//EquationData::param_omega*EquationData::param_sigma(i)*EquationData::constant_mu0;
        EquationData::param_kappa_im(i) = EquationData::param_sigma(i)*EquationData::param_omega*EquationData::constant_mu0;
      }
      else
      {
        EquationData::param_kappa_re(i) = EquationData::param_regularisation*EquationData::param_omega*EquationData::constant_mu0;
        EquationData::param_kappa_im(i) = 0.0;//EquationData::param_regularisation*EquationData::param_omega*EquationData::constant_mu0;
        //0.0;//EquationData::param_regularisation*EquationData::constant_mu0*EquationData::param_omega;
      }
    }    
  }
  template <int dim>
  void cubeBenchmark<dim>::process_mesh(bool neumann_flag)
  {
    // Routine to process the read in mesh
    // Here we set the boundary ids for the boundary conditions
    // and may choose to mark the boundary as curved, etc.
    // Can also perform mesh refinement here.
    
    // Set boundaries to neumann (boundary_id = 10)
    typename Triangulation<dim>::active_cell_iterator
      cell = tria.begin (),
      endc = tria.end();
    if (neumann_flag)
    {
      for (; cell!=endc; ++cell)
      {
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
          if (cell->face(face)->at_boundary())
          {
            cell->face(face)->set_boundary_indicator (10);
          }
        }
      }
    }
    
    // refine the inner mesh
//     cell = tria.begin_active();
//     for (; cell!=endc; ++cell)
//     {
//       if (cell->material_id() == 1)
//       {
//         cell->set_refine_flag();
//       }                         
//     }
//     tria.execute_coarsening_and_refinement ();

    // refine around inner mesh
//     double threshold = 0.05;
//     for (unsigned int r=0; r < 3; ++r)
//     {
//       cell = tria.begin_active();
//       for (; cell!=endc; ++cell)
//       {
//         Point<dim> cell_centre = cell->center();
//         if (sqrt(cell_centre.square()) < threshold )
//         {
//           cell->set_refine_flag();
//         }                         
//       }
//       tria.execute_coarsening_and_refinement ();    
//       threshold=threshold/2;
//     }
    

    // refine the outer mesh
/*    cell = tria.begin_active();
    for (; cell!=endc; ++cell)
    {
      if (cell->material_id() == 0)
      {
        cell->set_refine_flag();
      }                         
    }
    tria.execute_coarsening_and_refinement ();  */  

  }
  template <int dim>
  void cubeBenchmark<dim>::run(std::string input_filename, 
                                 std::string output_filename)
  {
    
    ParameterHandler prm;
    InputTools::ParameterReader param(prm);    
    param.read_and_copy_parameters(input_filename);
    
    InputTools::read_in_mesh<dim>(IO_Data::mesh_filename,
                                  tria);
    
    process_mesh(true);
    initialise_materials();
    
    
    std::cout << "Number of active cells:       "
    << tria.n_active_cells()
    << std::endl;
    
    // Now setup the forward problem:
    dof_handler.distribute_dofs (fe);
    ForwardSolver::EddyCurrent<dim, DoFHandler<dim>> eddy(dof_handler,
                                                          fe,
                                                          PreconditionerData::use_direct);
    
    std::cout << "Number of degrees of freedom: "
    << dof_handler.n_dofs()
    << std::endl;
    
    // assemble the matrix for the eddy current problem:
    std::cout << "Assembling System Matrix...." << std::endl;
    eddy.assemble_matrices(dof_handler);
    std::cout << "Matrix Assembly complete. " << std::endl;
    
    // initialise the linear solver - precomputes any inverses for the preconditioner, etc:
    std::cout << "Initialising Solver..." << std::endl;
    eddy.initialise_solver();
    std::cout << "Solver initialisation complete. " << std::endl;

    // Now solve for each excitation coil:
    std::cout << "Solving... " << std::endl;
    Vector<double> solution;
    
    // construct RHS for this field:
    Vector<double> uniform_field(dim+dim);
    uniform_field=0;
    uniform_field(2)=1;
    backgroundField::conductingCube<dim> boundary_conditions(uniform_field);
    
    // Check the exact field:
    /*
    std::vector<Vector<double> > sol(2, Vector<double>(dim+dim));
    std::vector<Point<dim> > points(2);
    points[0]=Point<dim> (0.6,0.6,0.6);
    points[1]=Point<dim> (0.25,0.25,0.25);
    boundary_conditions.vector_value_list(points, sol);
    for (unsigned int k=0; k<2; ++k)
    {
      for (unsigned int i=0; i<dim; ++i)
      {
        std::cout << sol[k](i) << " + " << sol[k](i+dim) << "*i" << std::endl;
      }
    }
    */
    
    // assemble rhs
    eddy.assemble_rhs(dof_handler,
                      boundary_conditions);
    
    // solve system & storage in the vector of solutions:
    eddy.solve(solution);

    std::cout << "Computed solution. " << std::endl;
    
//     {
//       std::ostringstream tmp;
//       tmp << output_filename;    
//       OutputTools::output_to_vtk<dim, DoFHandler<dim>>(dof_handler,
//                                                        solution,
//                                                        tmp.str(),
//                                                        boundary_conditions);
//     }
    
    // Output the solution fields along a line to a text file:
    Point <dim> xaxis;
    Point <dim> yaxis;
    Point <dim> zaxis;
    Point <dim> diagaxis;
    if (MeshData::boundary_shape=="cube")
    { 
      xaxis = Point<dim> (MeshData::xmax, 0.0, 0.0);
      yaxis = Point<dim> (0.0, MeshData::ymax, 0.0);
      zaxis = Point<dim> (0.0, 0.0, MeshData::zmax);
      diagaxis = Point<dim> (MeshData::xmax, MeshData::ymax, MeshData::zmax);
    }
    else if (MeshData::boundary_shape=="cylinder_z")
    {
      xaxis = Point<dim> (MeshData::radius, 0.0, 0.0);
      yaxis = Point<dim> (0.0, MeshData::radius, 0.0);
      zaxis = Point<dim> (0.0, 0.0, MeshData::zmax);
      diagaxis = Point<dim> (MeshData::radius*cos(numbers::PI/4.0), MeshData::radius*sin(numbers::PI/4.0), MeshData::zmax);
    }
    
    {
    std::stringstream tmp;
    tmp << output_filename << "_xaxis";
    std::string xaxis_str = tmp.str();
    OutputTools::output_radial_values<dim> (dof_handler,
                                       solution,
                                       boundary_conditions,
                                       uniform_field,
                                       xaxis,
                                       xaxis_str);
    }
    {
      std::stringstream tmp;
      tmp << output_filename << "_yaxis";
      std::string yaxis_str = tmp.str();
      OutputTools::output_radial_values<dim> (dof_handler,
                                              solution,
                                              boundary_conditions,
                                              uniform_field,
                                              yaxis,
                                              yaxis_str);
    }
    {
      std::stringstream tmp;
      tmp << output_filename << "_zaxis";
      std::string zaxis_str = tmp.str();
      OutputTools::output_radial_values<dim> (dof_handler,
                                       solution,
                                       boundary_conditions,
                                       uniform_field,
                                       zaxis,
                                       zaxis_str);
    }
    {
      std::stringstream tmp;
      tmp << output_filename << "_diagaxis";
      std::string diagaxis_str = tmp.str();
      OutputTools::output_radial_values<dim> (dof_handler,
                                              solution,
                                              boundary_conditions,
                                              uniform_field,
                                              diagaxis,
                                              diagaxis_str);
    }
    
    /*
    {
      // TODO: automate the start/end point process
      //        and make the direction an input variable.
      unsigned int n_points = 50;
      std::vector<Point<dim>> measurement_points(n_points);
      const double start_point=0.01/2+0.0001;
      const double end_point = 0.02;
      const double inc = (end_point-start_point)/n_points; // deliberately miss out boundary point as it's set by BCs
      for (unsigned int i=0; i<measurement_points.size(); ++i)
      {
        measurement_points[i] = Point<dim> (start_point+i*inc,
                                            0.0,//start_point+i*inc,
                                            0.0);//start_point+i*inc);
      }
      std::vector<Vector<double>> field_values_exact(measurement_points.size(),
                                                     Vector<double> (dim+dim));
      std::vector<Vector<double>> field_values_approx(measurement_points.size(),
                                                      Vector<double> (dim+dim));
      
      std::vector<Vector<double>> curl_values_exact(measurement_points.size(),
                                                     Vector<double> (dim+dim));
      std::vector<Vector<double>> curl_values_approx(measurement_points.size(),
                                                      Vector<double> (dim+dim));      
      
      std::vector<Vector<double>> perturbed_field_values_exact(measurement_points.size(),
                                                               Vector<double> (dim+dim));
      std::vector<Vector<double>> perturbed_field_values_approx(measurement_points.size(),
                                                                Vector<double> (dim+dim));
      
      boundary_conditions.vector_value_list(measurement_points,
                                            field_values_exact);                                      
      boundary_conditions.curl_value_list(measurement_points,
                                          curl_values_exact);
      boundary_conditions.perturbed_field_value_list(measurement_points,
                                                     perturbed_field_values_exact);
      
      for (unsigned int i=0; i<measurement_points.size(); ++i)
      {
        VectorTools::point_value(dof_handler,
                                 solution,
                                 measurement_points[i],
                                 field_values_approx[i]);
        
        Vector<double> temp_curl(dim+dim);
        MyVectorTools::point_curl(dof_handler,
                                  solution,
                                  measurement_points[i],
                                  temp_curl);
        for (unsigned int d=0; d<dim+dim; ++d)
        {
          curl_values_approx[i](d) = temp_curl(d);
          perturbed_field_values_approx[i](d) = temp_curl(d) - uniform_field(d);
        }
      }
      
      // output to file:
      std::ostringstream tmp;
      tmp << output_filename << "_ptfield_p" << p_order << ".out";
      std::ofstream file(tmp.str());
      file.precision(32);
      for (unsigned int i=0; i<measurement_points.size(); ++i)
      {
        double r = sqrt(measurement_points[i].square());
        file << r << " ";
        // approx field:
        for (unsigned int d=0;d<dim+dim; ++d)
        {
          file << field_values_approx[i](d) << " ";
        }
        // exact field:
        for (unsigned int d=0;d<dim+dim; ++d)
        {
          file << field_values_exact[i](d) << " ";
        }
        
        // approx curl
        for (unsigned int d=0;d<dim+dim; ++d)
        {
          file << curl_values_approx[i](d) << " ";
        }
        // exact curl
        for (unsigned int d=0;d<dim+dim; ++d)
        {
          file << curl_values_exact[i](d) << " ";
        }
        
        // approx curl perturbed
        for (unsigned int d=0;d<dim+dim; ++d)
        {
          file << perturbed_field_values_approx[i](d) << " ";
        }
        // exact curl perturbed
        for (unsigned int d=0;d<dim+dim; ++d)
        {
          file << perturbed_field_values_exact[i](d) << " ";
        }        

        file << std::endl;
      }
      file.close();
    }*/
  }
}

int main (int argc, char* argv[])
{
  using namespace dealii;
  
  unsigned int dim = 3;
  // Set default input:
  unsigned int p_order = 0;
  std::string output_filename = "cube";
  std::string input_filename = "../input_files/sphere_benchmark.prm";
  
  // Allow for input from command line:
  if (argc > 0)
  {
    for (int i=1;i<argc;i++)
    {
      if (i+1 != argc)
      {
        std::string input = argv[i];
        if (input == "-p")
        {
          std::stringstream strValue;
          strValue << argv[i+1];
          strValue >> p_order;
        }
        if (input == "-i")
        {
          input_filename = argv[i+1];
        }
        if (input == "-o")
        {
          output_filename = argv[i+1];
        }
      }
    }
  }
  
  // Only output to logfile, not console:
  deallog.depth_console(0);
  std::ostringstream deallog_filename;
  deallog_filename << output_filename << "_p" << p_order << ".deallog";
  std::ofstream deallog_file(deallog_filename.str());
  deallog.attach(deallog_file);
  
  cubeBenchmark::cubeBenchmark<3> eddy_voltages(p_order);
  eddy_voltages.run(input_filename,
                    output_filename);
  
  deallog_file.close();
  return 0;
}