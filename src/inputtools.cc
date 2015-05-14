#include <inputtools.h>

using namespace dealii;

namespace InputTools
{
  template <int dim>
  void read_in_mesh (std::string mesh_name,
                     Triangulation<dim> &triangulation)
                                              
  {
    // Intended for reading in .ucd files generated by Cubit/Trelis
    // These meshes have material ids (blocks) starting at 1
    // so we take 1 away from all material_id's in the mesh.
    
    GridIn<dim> gridin;
    gridin.attach_triangulation(triangulation);
    std::ifstream mesh_file(mesh_name);
    gridin.read_ucd(mesh_file);

    // Adjust material_id to start at 0 instead of 1.
    for (typename Triangulation<dim>::active_cell_iterator
      cell = triangulation.begin_active();
    cell != triangulation.end();
    ++cell)
    {
      cell->set_material_id(cell->material_id()-1);
    }  
  }
  // Template instantiation
  template void read_in_mesh<3>(std::string,
                                Triangulation<3> &);
  
  
  // CLASS PARAMETERREADER MEMBERS
  ParameterReader::ParameterReader(ParameterHandler &paramhandler)
  :
  prm(paramhandler)
  {}
  
  ParameterReader::~ParameterReader ()
  {}
  
  void ParameterReader::declare_parameters()
  {
    // Declare subsections with:
    //      prm.enter_subsection("subsection_name")
    // declare entries within subsection with:
    //       prm.declare("entry_name", "default_value", Pattern::type, "Description")
    
    prm.enter_subsection ("Mesh Data");
    {
      prm.declare_entry("external mesh", "true",
                        Patterns::Bool(),
                        "Use an external mesh file (specified by mesh file).");      
      
      prm.declare_entry("mesh file", "mesh.ucd",
                        Patterns::Anything(),
                        "Name of the mesh file (with extension)");
      
      prm.declare_entry("recovery region id", "1",
                        Patterns::Integer(0,1000),
                        "Recovery region material ID");
      
      prm.declare_entry("background region id", "0",
                        Patterns::Integer(0,1000),
                        "Background region material ID (no recovery performed here)");
      
      // Geometry info about mesh:
      // Can pass info to the code so we can use particular shapes & geometry
      // information below:
      prm.declare_entry("boundary shape", "unspecified",
                        Patterns::Anything(),
                        "Shape of boundary (cube/cylinder_z)");
      
      // general & box shapes
      // Typically assume centred at (0,0,0), but allow it to change.
      prm.declare_entry("xmax", "0.0",
                        Patterns::Double(),
                        "Maximum value of x");
      
      prm.declare_entry("xmin", "0.0",
                        Patterns::Double(),
                        "Minimum value of x");
      
      prm.declare_entry("ymax", "0.0",
                        Patterns::Double(),
                        "Maximum value of y");
      
      prm.declare_entry("ymin", "0.0",
                        Patterns::Double(),
                        "Minimum value of y");
      
      prm.declare_entry("zmax", "0.0",
                        Patterns::Double(),
                        "Maximum value of z");
      
      prm.declare_entry("zmin", "0.0",
                        Patterns::Double(),
                        "Minimum value of z");
      
      prm.declare_entry("centre","0.0, 0.0, 0.0",
                        Patterns::List(Patterns::Double(),3,3,","),
                        "Centre of mesh (if applicable)");
      
      // Cylinder shapes:
      // Can also use zmax for cylinder
      prm.declare_entry("radius", "0.0",
                        Patterns::Double(0),
                        "Radius of mesh (if applicable)");
      
      prm.declare_entry("height", "0.0",
                        Patterns::Double(0),
                        "Height of mesh (if applicable)");      
      
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("Output Parameters");
    {
      prm.declare_entry("Output filename", "solution",
                        Patterns::Anything(),
                        "Name of the output file (without extension)");
      
      prm.declare_entry("Output filetype", "vtk",
                        Patterns::Anything(),
                        "Output file extension");
      
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("Material Parameters");
    {
      prm.declare_entry("omega", "1.0",
                        Patterns::Double(0),
                        "Angular frequency");
      
      prm.declare_entry("regularisation parameter", "0.01",
                        Patterns::Double(0),
                        "Regularisation Parameter");
      
      prm.declare_entry("background epsilon", "0.0",
                        Patterns::Double(0),
                        "Background permittivity");
      
      prm.declare_entry("background mur", "1.0",
                        Patterns::Double(0),
                        "Background (relative) permeability");
      
      prm.declare_entry("background sigma", "0.0",
                        Patterns::Double(0),
                        "Background conductivity");
      
      prm.declare_entry("object epsilon", "0.0",
                        Patterns::Double(0),
                        "Object permittivity");
      
      prm.declare_entry("object mur", "1.0",
                        Patterns::Double(0),
                        "Object (relative) permability");
      
      prm.declare_entry("object sigma", "0.0",
                        Patterns::Double(0),
                        "Object conductivity");
      
      
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("Excitation Coil Data");
    {
      // TODO:
      // For now we can only handle an array of excitation coils
      // which lie in a circle around a particular point. These are
      // then used to generate a set of positions and directions.
      // Ideally we would be able to read in a list of
      // coil positions and directions. This would mean we would
      // only read in the number of coils and their position/directions.
      
      prm.declare_entry("number of coils", "1",
                        Patterns::Integer(1,100),
                        "Number of excitation coils in the array");
      
      prm.declare_entry("array centre","0.0, 0.0, 0.0",
                        Patterns::List(Patterns::Double(),3,3,","),
                        "Centre of the excitation coil array");
      
      prm.declare_entry("array radius", "0.0145",
                        Patterns::Double(0),
                        "Radius of the excitation coil array"); 
                       
      prm.declare_entry("array start angle", "0.0",
                        Patterns::Double(),
                        "Angle (in degrees) to the position of the first excitation coil");
      
      prm.declare_entry("coil radius", "0.025",
                       Patterns::Double(0),
                       "Radius of an excitation coil");                       
    }
    prm.leave_subsection ();
    
    prm.enter_subsection ("Sensor Coil Data");
    {
      // TODO:
      // For now we can only handle an array of sensor coils
      // which lie in a circle around a particular point. These are
      // then used to generate a set of positions and directions.
      // Ideally we would be able to read in a list of
      // coil positions and directions. This would mean we would
      // only read in the number of coils and their position/directions.
      // Not sure how easy this is to do when the list length is unknown??
      
      prm.declare_entry("number of coils", "1",
                        Patterns::Integer(1,100),
                        "Number of sensor coils");
      
      prm.declare_entry("array centre","0.0, 0.0, 0.0",
                        Patterns::List(Patterns::Double(),3,3,","),
                        "Centre of the sensor coil array");
     
      prm.declare_entry("array radius", "0.0135",
                        Patterns::Double(0),
                        "Radius of the sensor coil array");
      
      prm.declare_entry("array start angle", "0.0",
                        Patterns::Double(),
                        "Angle (in degrees) to the position of the first sensor coil");
      
      prm.declare_entry("coil radius", "0.025",
                       Patterns::Double(0),
                       "Radius of an individual sensor coil");
    }
    prm.leave_subsection ();
    
    prm.enter_subsection("Preconditioner Data");
    {
      // Disable the gmres solver:
      prm.declare_entry("use direct", "false",
                        Patterns::Bool(),
                        "Enable the sparse direct solver (disables the GMRES solver)");
      // Data for the preconditioner
      prm.declare_entry("diagonal strengthening", "1.0",
                        Patterns::Double(0),
                        "SparseILU diagonal strengthening");
      
      prm.declare_entry("extra off diagonals", "0",
                        Patterns::Integer(0),
                        "SparseILU extra off diagonals");
    }
    prm.leave_subsection ();
    
    prm.enter_subsection("Polarization Tensor");
    {
      /* Read in via list: */
      prm.declare_entry("Polarization Tensor Real",
                        "1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0",
                        Patterns::List(Patterns::List(Patterns::Double(),3,3,","),3,3,";"),
                        "Real part of Polarization Tensor");
            
      prm.declare_entry("Polarization Tensor Imaginary",
                        "1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0",
                        Patterns::List(Patterns::List(Patterns::Double(),3,3,","),3,3,";"),
                        "Imaginary part of Polarization Tensor");
    }
    prm.leave_subsection ();
    

  }
  
  void ParameterReader::get_matrix_from_list(std::string entry, FullMatrix<double> &matrix_out, unsigned int matrix_size)
  {
    /* Outputs a matrix read in by ParameterHandler as a List (Patterns::List(Patterns::Double())
     * - Assumes a square matrix but could be extended to handle non-square.
     * - Also assumes the separator is a comma - could be extended to handle any.
     * 
     * Requires that the matrix is entered using ; to signal a new row
     * and , to signal a new column
     */
    
    std::stringstream wholeline(prm.get(entry));
    std::string rowStr;
    std::string val;
    for (unsigned int i=0;i<matrix_size;i++)
    {
      std::getline(wholeline,rowStr,';');
      std::stringstream row(rowStr);
      for (unsigned int j=0;j<matrix_size;j++)
      {
        std::getline(row, val, ',');
        std::stringstream converter(val);
        converter >> matrix_out(i,j);
      }
    }
  }
  
  void ParameterReader::get_vector_from_list(std::string entry, Vector<double> &vector_out, unsigned int vector_length)
  {
    /* Outputs a vector read in by ParameterHandler as a List (Patterns::List(Patterns::Double())
     */
    
    std::stringstream wholeline(prm.get(entry));
    std::string val;
    for (unsigned int i=0;i<vector_length;i++)
    {
      std::getline(wholeline,val,',');
      std::stringstream converter(val);
      converter >> vector_out(i);
    }
  }
  
  
  void ParameterReader::copy_to_equation_data()
  {
    // temp vector data:
    Vector<double> temp_vector1(3);
    Vector<double> temp_vector2(3);
    
    // Input data:
    prm.enter_subsection("Mesh Data");
    
    MeshData::external_mesh = prm.get_bool("external mesh");
    //TODO: remove and leave only MeshData::mesh_filename (will have to update all codes first.)
    IO_Data::mesh_filename = prm.get("mesh file");
    MeshData::mesh_filename = prm.get("mesh file");
    
    InverseProblemData::recovery_region_id = prm.get_integer("recovery region id");
    InverseProblemData::background_region_id = prm.get_integer("background region id");
    
    // Geometry info about mesh:
    MeshData::boundary_shape = prm.get("boundary shape");
    MeshData::xmax = prm.get_double("xmax");
    MeshData::xmin = prm.get_double("xmin");

    MeshData::ymax = prm.get_double("ymax");
    MeshData::ymin = prm.get_double("ymin");
      
    MeshData::zmax = prm.get_double("zmax");
    MeshData::zmin = prm.get_double("zmin");
      
    MeshData::height = prm.get_double("radius");
    MeshData::radius = prm.get_double("height");
    get_vector_from_list("centre", temp_vector1, 3);
    for (unsigned int d=0; d<3; ++d)
    {
      MeshData::centre(d) = temp_vector1(d);
    }
    prm.leave_subsection();
    
    // Output data:
    prm.enter_subsection("Output Parameters");
    
    IO_Data::output_filename = prm.get("Output filename");
    IO_Data::output_filetype = prm.get("Output filetype");
    
    prm.leave_subsection();
    
    // Material parameters:
    prm.enter_subsection("Material Parameters");
    
    EquationData::param_epsilon_background = prm.get_double("background epsilon");
    EquationData::param_mur_background = prm.get_double("background mur");
    EquationData::param_sigma_background = prm.get_double("background sigma");
    
    EquationData::param_omega = prm.get_double("omega");
    
    EquationData::param_regularisation = prm.get_double("regularisation parameter");
    
    EquationData::param_epsilon_conducting = prm.get_double("object epsilon");
    EquationData::param_mur_conducting = prm.get_double("object mur");
    EquationData::param_sigma_conducting = prm.get_double("object sigma");
    
    prm.leave_subsection();
    
    // Excitation and Sensor Coil Data:
    // Data:
    //      number of coils
    //      array centre
    //      array radius
    //      array start angle
    //      coil radius
    
    prm.enter_subsection("Excitation Coil Data");
    
    ExcitationCoilData::number_of_coils = prm.get_integer("number of coils");
    
    get_vector_from_list("array centre", temp_vector1, 3);
    for (unsigned int i=0; i<3; ++i)
    {
      ExcitationCoilData::array_centre[i] = temp_vector1(i);
    }
    ExcitationCoilData::array_radius = prm.get_double("array radius");
    ExcitationCoilData::array_angle = prm.get_double("array start angle");
    ExcitationCoilData::coil_radius = prm.get_double("coil radius");
    
    prm.leave_subsection();
    
    prm.enter_subsection("Sensor Coil Data");
    
    SensorCoilData::number_of_coils = prm.get_integer("number of coils");
    
    get_vector_from_list("array centre", temp_vector1, 3);
    for (unsigned int i=0; i<3; ++i)
    {
      SensorCoilData::array_centre[i] = temp_vector1(i);
    }
    SensorCoilData::array_radius = prm.get_double("array radius");
    SensorCoilData::array_angle = prm.get_double("array start angle");
    SensorCoilData::coil_radius = prm.get_double("coil radius");
    
    prm.leave_subsection();
    
    prm.enter_subsection("Preconditioner Data");
    PreconditionerData::use_direct = prm.get_bool("use direct");
    PreconditionerData::strengthen_diagonal = prm.get_double("diagonal strengthening");
    PreconditionerData::extra_off_diagonals = prm.get_integer("extra off diagonals");
    prm.leave_subsection();
    
    prm.enter_subsection("Polarization Tensor");
    get_matrix_from_list("Polarization Tensor Real", PolarizationTensor::polarizationTensor_re,3);
    get_matrix_from_list("Polarization Tensor Imaginary", PolarizationTensor::polarizationTensor_im,3);
    prm.leave_subsection();
  }
  
  void ParameterReader::read_parameters (const std::string parameter_file)
  {
    declare_parameters();
    prm.read_input (parameter_file);
  }
  void ParameterReader::read_and_copy_parameters (const std::string parameter_file)
  {
    declare_parameters();
    prm.read_input (parameter_file);
    copy_to_equation_data();
  }
  void ParameterReader::copy_parameters ()
  {
    copy_to_equation_data();
  }
  // template instantiation not needed (derived class) ???
  // END CLASS PARAMETERREADER
}