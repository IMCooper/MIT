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

namespace wavepropBenchmark
{
  
  template <int dim>
  class wavepropBenchmark
  {
  public:
    wavepropBenchmark (const unsigned int order);
    ~wavepropBenchmark ();
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
  wavepropBenchmark<dim>::wavepropBenchmark(const unsigned int order)
  :
  fe (FE_Nedelec<dim>(order), 2),
  dof_handler (tria),
  p_order(order)
  {
  }
  
  template <int dim>
  wavepropBenchmark<dim>::~wavepropBenchmark ()
  {
    dof_handler.clear ();  
  }
  
  template <int dim>
  void wavepropBenchmark<dim>::initialise_materials()
  {
    // Allow for testing of eddy current input files:
    // Essentially, we discount the sigma/eps values
    // but the mur values are used.
    EquationData::param_mur.reinit(2);
    EquationData::param_mur(0) = EquationData::param_mur_conducting;       
    EquationData::param_mur(1) = EquationData::param_mur_conducting;
    
    EquationData::param_sigma.reinit(2);
    EquationData::param_sigma(0) = EquationData::param_sigma_conducting;
    EquationData::param_sigma(1) = EquationData::param_sigma_conducting;
    
    EquationData::param_epsilon.reinit(2);
    EquationData::param_epsilon(0) = EquationData::param_epsilon_conducting;
    EquationData::param_epsilon(1) = EquationData::param_epsilon_conducting;
    
    
    // kappa = -omega^2*epr + i*omega*sigma;
    // i.e. kappa_re = -omega^2*epsilon
    //      kappa_im = omega*sigma
    EquationData::param_kappa_re.reinit(EquationData::param_mur.size());
    EquationData::param_kappa_im.reinit(EquationData::param_mur.size());
    for (unsigned int i=0;i<EquationData::param_mur.size();i++) // note mur and kappa must have same size:
    {
      
      EquationData::param_kappa_re(i) = -EquationData::param_omega*EquationData::param_omega;
      EquationData::param_kappa_im(i) = 0.0;
      
    }    
  }
  template <int dim>
  void wavepropBenchmark<dim>::process_mesh(bool neumann_flag)
  {
    // Routine to process the read in mesh
    // Here we set the boundary ids for the boundary conditions
    // and may choose to mark the boundary as curved, etc.
    // Can also perform mesh refinement here.
    
    // Set boundaries to neumann (boundary_id = 10)
    if (neumann_flag)
    {
      typename Triangulation<dim>::active_cell_iterator
      cell = tria.begin (),
      endc = tria.end();
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
    

    // Now refine the outer mesh
    /*
    cell = tria.begin_active();
    for (; cell!=endc; ++cell)
    {
//       if (cell->material_id() == 0)
//       {
        cell->set_refine_flag();
//       }                         
    }
    tria.execute_coarsening_and_refinement ();    
    */
  }
  template <int dim>
  void wavepropBenchmark<dim>::run(std::string input_filename, 
                                 std::string output_filename)
  {
    
    ParameterHandler prm;
    InputTools::ParameterReader param(prm);
    param.read_and_copy_parameters(input_filename);
    

    /*
    InputTools::read_in_mesh<dim>(IO_Data::mesh_filename,
                                  tria);
                                  */

    
    GridGenerator::hyper_cube (tria, -0.0125, 0.0125);
    tria.refine_global (2);
    GridTools::distort_random (0.2, tria, false);
    
    process_mesh(false);
    
    // construct RHS for this field:
    Vector<double> p_wave(dim);
    Vector<double> k_wave(dim);
    // along x axis
    k_wave(0) = 0.0;
    k_wave(1) = 0.0;
    k_wave(2) = 1.0;
    p_wave(0) = 1.0;
    p_wave(1) = 0.0;
    p_wave(2) = 0.0;

    // diagonal:
//     k_wave(0) = -sqrt(3)/2;
//     k_wave(1) = -sqrt(3)/2;
//     k_wave(2) = sqrt(3);
//     p_wave(0) = 1./sqrt(3.);
//     p_wave(1) = 1./sqrt(3.);
//     p_wave(2) = 1./sqrt(3.);
    backgroundField::WavePropagation<dim> boundary_conditions(k_wave,
                                                              p_wave);
    
    // Override any other value of omega here
    // as it must be such that omega = |k|
    EquationData::param_omega = k_wave.l2_norm();
    
    
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
    /*
    {
      std::ostringstream tmp;
      tmp << output_filename;    
      OutputTools::output_to_vtk<dim, DoFHandler<dim>>(dof_handler,
                                                       solution,
                                                       tmp.str(),
                                                       boundary_conditions);
    }
    */
    
    // Output the perturbed field to a text file:
    
    {
      unsigned int n_points = 100;
      std::vector<Point<dim>> measurement_points(n_points);
      Point<dim> end_point (MeshData::xmax,MeshData::ymax,MeshData::zmax); // Eventually make more changeable.
      const double inc = 1.0/n_points;
      for (unsigned int i=0; i<measurement_points.size(); ++i)
      {
        measurement_points[i] = i*inc*end_point;
      }
      std::vector<Vector<double>> field_values_exact(measurement_points.size(),
                                                     Vector<double> (dim+dim));
      std::vector<Vector<double>> field_values_approx(measurement_points.size(),
                                                      Vector<double> (dim+dim));
      
      std::vector<Vector<double>> curl_values_exact(measurement_points.size(),
                                                    Vector<double> (dim+dim));
      std::vector<Vector<double>> curl_values_approx(measurement_points.size(),
                                                     Vector<double> (dim+dim));
      
      boundary_conditions.vector_value_list(measurement_points,
                                            field_values_exact);                                      
      boundary_conditions.curl_value_list(measurement_points,
                                          curl_values_exact);
      
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
        }
      }
      // output to file:
      std::ostringstream tmp;
      tmp << output_filename << "_fields_p" << p_order << ".out";
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
        
        file << std::endl;
      }
      file.close();
    }
    // Output error to screen:
    double hcurlerr = MyVectorTools::calcErrorHcurlNorm(dof_handler,
                                                        solution,
                                                        boundary_conditions);
    std::cout << "HCurl Error: " << hcurlerr << std::endl;
  }
}

int main (int argc, char* argv[])
{
//  using namespace dealii;
  
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
  
  wavepropBenchmark::wavepropBenchmark<3> eddy_voltages(p_order);
  eddy_voltages.run(input_filename,
                    output_filename);
  
  deallog_file.close();
  return 0;
}
