#include <backgroundfield.h>

using namespace dealii;

namespace backgroundField
{
  // DIPOLESOURCE
  template<int dim>
  DipoleSource<dim>::DipoleSource(const Point<dim> &input_source_point,
                                  const Tensor<1, dim> &input_coil_direction)
  :
  source_point(input_source_point),
  coil_direction(input_coil_direction)
  {}
  
  // members:
  template <int dim>
  void DipoleSource<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                             std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    const unsigned int n_points = points.size();       
        
    Tensor<1,dim> shifted_point;
    Tensor<1,dim> result; 
    for (unsigned int k=0; k<points.size(); ++k)
    {
      const Point<dim> &p = points[k];
      /* Work out the vector (stored as a tensor so we can use cross_product)
       * from the source point to the current point, p
       */
      for (unsigned int i = 0; i < dim; ++i)
      {
        shifted_point[i] = p(i) - source_point(i);
      }
      double rad = p.distance(source_point);
      double factor = EquationData::constant_mu0*1.0/(4.0*numbers::PI*rad*rad*rad);
      
      cross_product(result, coil_direction, shifted_point);
      result *= factor;
      for (unsigned int i = 0; i < dim; ++i)
      {
        // Real
        value_list[k](i) = result[i];
        // Imaginary
        value_list[k](i+dim) = 0.0;
      }
    }
  }
  template <int dim>
  void DipoleSource<dim>::curl_value_list (const std::vector<Point<dim> > &points,
                                           std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    const unsigned int n_points = points.size();
    
    FullMatrix<double> rhat(dim);
    FullMatrix<double> D2G(dim);
    FullMatrix<double> eye(IdentityMatrix(3));
    Vector<double> result(dim);
    // create a vector-version of the tensor coil_direction.
    // There may be a far better way to deal with this.... (TODO)
    // e.g. Use Tensor<2,dim> for the matrices instead, making the whole thing more tensor based.
    Vector<double> coil_direction(dim);
    for (unsigned int i = 0; i < dim; ++i)
    {
      coil_direction(i) = coil_direction[i];
    }

    Tensor<1,dim> shifted_point;
    Vector<double> scaled_vector(dim);

    for (unsigned int k=0; k<points.size(); ++k)
    {
      const Point<dim> &p = points[k];
      // Work out the vector (stored as a tensor so we can use cross_product)
      // from the source point to the current point, p
      for (unsigned int i = 0; i < dim; ++i)
      {
        shifted_point[i] = p(i) - source_point(i);
      }
      double rad = p.distance(source_point);
      double factor = 1.0/(4.0*numbers::PI*rad*rad*rad);

      // Construct D2G
      for (unsigned int i = 0; i < dim; ++i)
      {
        scaled_vector(i)=shifted_point[i]/rad;
      }
      rhat.outer_product(scaled_vector,scaled_vector);
      D2G=0;
      D2G.add(3.0,rhat,-1.0,eye);

      D2G.vmult(result, coil_direction);

      result *= factor;
      for (unsigned int i=0;i<dim;i++)
      {
        // Real
        value_list[k](i) = result[i];
        // Imaginary
        value_list[k](i+dim) = 0.0;
      }
    }
  }
  template class DipoleSource<3>;
  // END DIPOLESOURCE
  
  // CONDUCTINGSPHERE
  template<int dim>
  conductingSphere<dim>::conductingSphere(double sphere_radius,
                                          const Vector<double> &uniform_field)
  :
  sphere_radius(sphere_radius),
  uniform_field_re(dim),
  uniform_field_im(dim)
  {
    // copy over uniform field
    for (unsigned int d=0;d<dim;++d)
    {
      uniform_field_re(d)=uniform_field(d);
      uniform_field_im(d)=uniform_field(d+dim);
    }
    constant_B_magnitude = sqrt(uniform_field_re.norm_sqr() + uniform_field_im.norm_sqr());
    
    // calculate the required constants:
    // TODO
    double sigma = EquationData::param_sigma_conducting;
    double mu_c = EquationData::param_mur_conducting*EquationData::constant_mu0;
    double mu_n = EquationData::param_mur_background*EquationData::constant_mu0;
    double omega = EquationData::param_omega;
    constant_p = sigma*mu_c*omega;
    std::complex<double> temp(0, constant_p);
    std::complex<double> v = sqrt(temp)*sphere_radius;
    
    // Calculate some useful terms:
    temp = sqrt(2.0/(numbers::PI*v));
    std::complex<double> besselhalf_plus = temp*sinh(v);
    std::complex<double> besselhalf_minus = temp*cosh(v);
    
    // If conducting & non-conducting region have same mu, then it simplifies:
    
    if (EquationData::param_mur_background == EquationData::param_mur_conducting)
    {
      constant_C = 3.0*pow(sqrt(sphere_radius),3)/(v*besselhalf_plus);
      constant_D = pow(sphere_radius,3)*( 3.0*v*besselhalf_minus - (3.0+v*v)*besselhalf_plus )
                   / ( v*v*besselhalf_plus );
    }
    else
    {    
       constant_C = 3.0*mu_c*v*pow(sqrt(sphere_radius),3) / ( (mu_c-mu_n)*v*besselhalf_minus + (mu_n*(1.0+v*v)-mu_c)*besselhalf_plus );
       constant_D = pow(sphere_radius,3)*( (2.0*mu_c + mu_n)*v*besselhalf_minus - (mu_n*(1.0+v*v) + 2.0*mu_c)*besselhalf_plus )
                    / ( (mu_c - mu_n)*v*besselhalf_minus + (mu_n*(1.0+v*v)-mu_c)*besselhalf_plus );
    }
                    
    std::cout << constant_C.real() << " + "<<constant_C.imag() << "i" << std::endl;
    std::cout << constant_D.real() << " + "<<constant_D.imag() << "i" << std::endl;
  }
  
  template <int dim>
  void conductingSphere<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                             std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    const unsigned int n_points = points.size();       

    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      const Point<dim> &p = points[k];
      double r, theta, phi;
      std::complex<double> factor;
      // Convert (x,y,z) to (r,theta,phi), spherical polars.
      r = sqrt(p.square());
      theta = acos(p(2)/r);
//       theta = atan2(p(1),hypot(p(0),p(2)));
      phi = atan2(p(1),p(0)); // TODO
      
      if (r < 1e-4)
      {
        factor=0.0;
        value_list[k]=0;
      }        
      else if (r < sphere_radius)
      {
        std::complex<double> temp(0, constant_p);
        std::complex<double> v=sqrt(temp)*r;
        
        std::complex<double> bessel3halfs_plus = sqrt(2.0/(numbers::PI*v))*(cosh(v) - (1.0/v)*sinh(v));
        factor = 0.5*constant_B_magnitude*constant_C*bessel3halfs_plus*sin(theta)/sqrt(r);
      }
      else
      {
        factor = 0.5*constant_B_magnitude*(r + constant_D/(r*r))*sin(theta);
      }
      // Convert back to cartesian 
      // & split real/imaginary parts:
      value_list[k](0) = -factor.real()*sin(phi);
      value_list[k](1) = factor.real()*cos(phi);
      value_list[k](2) = 0.0;
      value_list[k](3) = -factor.imag()*sin(phi);
      value_list[k](4) = factor.imag()*cos(phi);
      value_list[k](5) = 0.0;
    }
  }
  
  template <int dim>
  void conductingSphere<dim>::curl_value_list (const std::vector<Point<dim> > &points,
                                             std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    
    // Use perturbed values:
    /*
    std::vector<Vector<double> > perturbed_values (value_list.size(), Vector<double> (dim+dim));    
    perturbed_field_value_list(points, perturbed_values);
    */
    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      
      // Use the perturbed values:
      /*
      for (unsigned int d=0; d<dim; ++d)
      {
        value_list[k](d) = perturbed_values[k](d) + uniform_field_re(d);
        value_list[k](d+dim) = perturbed_values[k](d+dim) + uniform_field_im(d);
      }
      */
      const Point<dim> &p = points[k];
      double r, theta, phi;
      // Convert (x,y,z) to (r,theta,phi), spherical polars.
      r = sqrt(p.square());
      theta = acos(p(2)/r);
//       theta = atan2(p(1),hypot(p(0),p(2)));
      phi = atan2(p(1),p(0)); // TODO
      
      if (r < 1e-4)
      {
        value_list[k]=0;
      }        
      else if (r < sphere_radius)
      {
        // TODO: confirm this is correct
        std::complex<double> temp(0, constant_p);
        std::complex<double> sqrt_ip = sqrt(temp);
        std::complex<double> v=sqrt(temp)*r;
        std::complex<double> bessel3halfs_plus = sqrt(2.0/(numbers::PI*v))*(cosh(v) - (1.0/v)*sinh(v));
        std::complex<double> besselhalf_plus = temp*sinh(v);
        
        std::complex<double> factor_r
        = (1/sqrt(r*r*r))*constant_B_magnitude*constant_C*bessel3halfs_plus*cos(theta);
        std::complex<double> factor_theta
        = -constant_B_magnitude*constant_C*sin(theta)/(4.0*r)
        *( (1/sqrt(r) + sqrt_ip )*bessel3halfs_plus + sqrt_ip*besselhalf_plus );
          
          // Convert to cartesian:
        std::complex<double> factor_x = factor_r*sin(theta)*cos(phi) + factor_theta*cos(theta)*cos(phi);
        std::complex<double> factor_y = factor_r*sin(theta)*sin(phi) + factor_theta*cos(theta)*sin(phi);
        std::complex<double> factor_z = factor_r*cos(theta) - factor_theta*sin(theta);
        
        value_list[k](0) = factor_x.real();
        value_list[k](1) = factor_y.real();
        value_list[k](2) = factor_z.real();
        value_list[k](3) = factor_x.imag();
        value_list[k](4) = factor_y.imag();
        value_list[k](5) = factor_z.imag();
      }
      else // r > sphere_radius
      {
        std::complex<double> factor_r = constant_B_magnitude*(1.0 + constant_D/(r*r*r))*cos(theta);
        std::complex<double> factor_theta = -constant_B_magnitude*(1.0 - constant_D/(2.0*r*r*r))*sin(theta);
        // Convert to cartesian:
        std::complex<double> factor_x = factor_r*sin(theta)*cos(phi) + factor_theta*cos(theta)*cos(phi);
        std::complex<double> factor_y = factor_r*sin(theta)*sin(phi) + factor_theta*cos(theta)*sin(phi);
        std::complex<double> factor_z = factor_r*cos(theta) - factor_theta*sin(theta);
        
        value_list[k](0) = factor_x.real();
        value_list[k](1) = factor_y.real();
        value_list[k](2) = factor_z.real();
        value_list[k](3) = factor_x.imag();
        value_list[k](4) = factor_y.imag();
        value_list[k](5) = factor_z.imag();
      }
    }
  }
  template<int dim>
  void conductingSphere<dim>::perturbed_field_value_list (const std::vector< Point<dim> > &points,
                                                          std::vector<Vector<double> > &value_list) const
  {
    // Returns the value of the perturbed field:
    // H_{p} = H - H_{0}
    //
    // TODO: Assume that the centre of the object is (0,0,0) for now
    
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));

    FullMatrix<double> rhat(dim);
    FullMatrix<double> D2G(dim);
    FullMatrix<double> eye(IdentityMatrix(3));
    Vector<double> result_re(dim);
    Vector<double> result_im(dim);
    
    Vector<double> scaled_vector(dim);
    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      const Point<dim> &p = points[k];

      double r = sqrt(p.square());
      if (r - sphere_radius < 0)
      {
        for (unsigned int i=0; i<dim; ++i)
        {
          // TODO: add the formula for the inside of the sphere.
          value_list[k](i) = 0.0;
          value_list[k](i+dim) = 0.0;
        }
      }
      else
      {        
        double factor = 1.0/(4.0*numbers::PI*r*r*r);
        // Construct D2G
        for (unsigned int i=0; i<dim; ++i)
        {
          scaled_vector(i)=p(i)/r;
        }
        rhat.outer_product(scaled_vector,scaled_vector);
        D2G=0;
        D2G.add(3.0,rhat,-1.0,eye);
        
        D2G.vmult(result_re, uniform_field_re);
        D2G.vmult(result_im, uniform_field_im); // D2G is real valued so no extra terms.
        
        result_re *= factor;
        result_im *= factor;
        // Now multiply by 2*pi*polarization_tensor
        // NOTE: the polarization tensor is diagonal: M = constant_D*identityMatrix.
        //       and it is complex valued.
        for (unsigned int i=0; i<dim; ++i)
        {
          value_list[k](i) = 2.0*numbers::PI*(constant_D.real()*result_re(i) - constant_D.imag()*result_im(i));
          value_list[k](i+dim) = 2.0*numbers::PI*(constant_D.imag()*result_re(i) + constant_D.real()*result_im(i));
        }
      }
    }
  }
  template<int dim>
  void conductingSphere<dim>::check_spherical_coordinates(const std::vector< Point<dim> > &points,
                                                          std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      const Point<dim> &p = points[k];
      double r, theta, phi;
      // Convert (x,y,z) to (r,theta,phi), spherical polars.
      r = sqrt(p.square());
//       theta = atan2(p(1),hypot(p(0),p(2)));
      theta = acos(p(2)/r);
      phi = atan2(p(1),p(0));
      value_list[k](0)=r;
      value_list[k](1)=theta;
      value_list[k](2)=phi;
    }
  }

  template class conductingSphere<3>; 
  // END CONDUCTINGSPHERE
  
  // CONDUCTINGCUBE
  template<int dim>
  conductingCube<dim>::conductingCube(const Vector<double> &uniform_field)
  :
  uniform_field_re(dim),
  uniform_field_im(dim)
  {
    // copy over uniform field
    for (unsigned int d=0;d<dim;++d)
    {
      uniform_field_re(d)=uniform_field(d);
      uniform_field_im(d)=uniform_field(d+dim);
    }
  }
  
  template <int dim>
  void conductingCube<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                               std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    const unsigned int n_points = points.size();       

    // No analytical solution, just set uniform far field conditions.
    for (unsigned int k=0; k<points.size(); ++k)
    {
      value_list[k]=0.0;
    }
  }
  
  template <int dim>
  void conductingCube<dim>::curl_value_list (const std::vector<Point<dim> > &points,
                                             std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    
    // No analytical solution, just set uniform far field conditions. 
    // TODO: switch over to using the perturbation tensor to calculate hte
    //       perturbed field, then add the uniform field.
    
    std::vector<Vector<double> > perturbed_value_list(value_list.size(), Vector<double> (dim+dim));
    perturbed_field_value_list(points, perturbed_value_list);
    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      for (unsigned int d=0; d<dim; ++d)
      {
        // Seems to only work when we treat the BC as n x curl E = 
        // TODO: work out why using curl(A) = B (computed from perturbation tensor formula) doesn't work.
        
        value_list[k](d) = perturbed_value_list[k](d) + uniform_field_re(d);
        //EquationData::constant_mu0*
        value_list[k](d+dim) = perturbed_value_list[k](d+dim) + uniform_field_im(d);
      }
    }
  }
  template<int dim>
  void conductingCube<dim>::perturbed_field_value_list (const std::vector< Point<dim> > &points,
                                                        std::vector<Vector<double> > &value_list) const
  {
    // Returns the value of the perturbed field:
    // H_{p} = H - H_{0}
    //
    // TODO: Assume that the centre of the object is (0,0,0) for now
    
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));

    FullMatrix<double> rhat(dim);
    FullMatrix<double> D2G(dim);
    FullMatrix<double> eye(IdentityMatrix(3));
    Vector<double> D2Gresult_re(dim);
    Vector<double> D2Gresult_im(dim);
    Vector<double> PTresult_re(dim);
    Vector<double> PTresult_im(dim);
    
    Vector<double> scaled_vector(dim);
    
    for (unsigned int k=0; k<points.size(); ++k)
    {
      const Point<dim> &p = points[k];

      double r = sqrt(p.square());
      if (r - cube_size < 0)
      {
        for (unsigned int i=0; i<dim; ++i)
        {
          value_list[k](i) = 0.0;
          value_list[k](i+dim) = 0.0;
        }
      }
      else
      {        
        double factor = 1.0/(4.0*numbers::PI*r*r*r);
        // Construct D2G
        for (unsigned int i=0; i<dim; ++i)
        {
          scaled_vector(i)=p(i)/r;
        }
        rhat.outer_product(scaled_vector,scaled_vector);
        D2G=0;
        D2G.add(3.0,rhat,-1.0,eye);
        
        D2G.vmult(D2Gresult_re, uniform_field_re);
        D2G.vmult(D2Gresult_im, uniform_field_im); // D2G is real valued so no extra terms.
        
        D2Gresult_re *= factor;
        D2Gresult_im *= factor;
        
        PolarizationTensor::polarizationTensor_re.vmult(PTresult_re, D2Gresult_re);
        PolarizationTensor::polarizationTensor_re.vmult(PTresult_im, D2Gresult_im);
        
        // Need to do imag*imag -> mult by -1
        D2Gresult_im *= -1.0;
        PolarizationTensor::polarizationTensor_im.vmult_add(PTresult_re, D2Gresult_im);
        PolarizationTensor::polarizationTensor_im.vmult_add(PTresult_im, D2Gresult_re);

        for (unsigned int i=0; i<dim; ++i)
        {
          value_list[k](i) = PTresult_re(i);
          value_list[k](i+dim) = PTresult_im(i);
        }
      }
    }
  }

  template class conductingCube<3>; 
  // END CONDUCTINGCUBE
  
  // WAVE PROPAGATION
  template<int dim>
  WavePropagation<dim>::WavePropagation(Vector<double> &k_wave,
                                        Vector<double> &p_wave)
  :
  k_wave(k_wave),
  p_wave(p_wave)
  {
    // Note: need k orthogonal to p
    //       with |k| < 1, |p| = 1.
    // Examples:
    // horizontal wave.
    //   k_wave(0) = 0.0;
    //   k_wave(1) = 0.1;
    //   k_wave(2) = 0.0;
    //   p_wave(0) = 1.0;
    //   p_wave(1) = 0.0;
    //   p_wave(2) = 0.0;
    // diagonal wave.
    //   k_wave(0) = -0.1;
    //   k_wave(1) = -0.1;
    //   k_wave(2) = 0.2;
    //   p_wave(0) = 1./sqrt(3.);
    //   p_wave(1) = 1./sqrt(3.);
    //   p_wave(2) = 1./sqrt(3.);
  
    // Make sure input is sane:
    Assert (k_wave.size() == 3, ExcDimensionMismatch (k_wave.size(), 3));
    Assert (k_wave.size() == k_wave.size(), ExcDimensionMismatch (k_wave.size(), k_wave.size()));
    const double delta = 1e-10;
    Assert (abs(p_wave.norm_sqr()-1.0) < delta, ExcIndexRange(p_wave.norm_sqr(), 1, 1));
  } 

  template <int dim>
  void WavePropagation<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                                std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    const unsigned int n_points = points.size();

    double exponent;
    for (unsigned int i=0; i<n_points; ++i)
    {
      const Point<dim> &p = points[i];
      exponent = k_wave(0)*p(0)+k_wave(1)*p(1)+k_wave(2)*p(2);
      for (unsigned int d=0; d<dim; ++d)
      {
        // Real:
        value_list[i](d) = ( p_wave(d) )*std::cos(exponent);
        // Imaginary:
        value_list[i](d+dim) = ( p_wave(d) )*std::sin(exponent);
      }
    }
  }
  template <int dim>
  void WavePropagation<dim>::curl_value_list(const std::vector<Point<dim> > &points,
                                             std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    const unsigned int n_points = points.size();
    
    
    double exponent;
    for (unsigned int i=0; i<n_points; ++i)
    {
      const Point<dim> &p = points[i];
      exponent = k_wave(0)*p(0)+k_wave(1)*p(1)+k_wave(2)*p(2);
      // Real:
      value_list[i](0) = -( k_wave(1)*p_wave(2) - k_wave(2)*p_wave(1) )*std::sin(exponent);
      value_list[i](1) = -( k_wave(2)*p_wave(0) - k_wave(0)*p_wave(2) )*std::sin(exponent);
      value_list[i](2) = -( k_wave(0)*p_wave(1) - k_wave(1)*p_wave(0) )*std::sin(exponent);
      // Imaginary:
      value_list[i](3) =  ( k_wave(1)*p_wave(2) - k_wave(2)*p_wave(1) )*std::cos(exponent);
      value_list[i](4) =  ( k_wave(2)*p_wave(0) - k_wave(0)*p_wave(2) )*std::cos(exponent);
      value_list[i](5) =  ( k_wave(0)*p_wave(1) - k_wave(1)*p_wave(0) )*std::cos(exponent);
    }
  }
  template class WavePropagation<3>; 
  // END WAVE PROPAGATION
  
  // POLYNOMIALTEST:
  template<int dim>
  polynomialTest<dim>::polynomialTest()//(Vector<double> &k_wave,
                                      //Vector<double> &p_wave)
  :
//   curlFunction<dim> (dim+dim)
  curlFunction<dim>()//dim+dim)
  // Don't know why I have to put these here, but won't compile without.. says
  // it's not of type curlFunction ???
//   k_wave(k_wave),
//   p_wave(p_wave)
  {}
  
  template <int dim>
  void polynomialTest<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                              std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    const unsigned int n_points = points.size();

    for (unsigned int i=0; i<n_points; ++i)
      {
        const Point<dim> &p = points[i];

        /* quadratic: */
        value_list[i](0) = p(0)*p(0);
        value_list[i](1) = p(1)*p(1);
        value_list[i](2) = p(2)*p(2);
        value_list[i](3) = p(0)*p(0);
        value_list[i](4) = p(1)*p(1);
        value_list[i](5) = p(2)*p(2);
      }
  }
  template <int dim>
  void polynomialTest<dim>::curl_value_list(const std::vector<Point<dim> > &points,
                                            std::vector<Vector<double> > &value_list) const
  {
    Assert(value_list.size() == points.size(), ExcDimensionMismatch(value_list.size(), points.size()));
    const unsigned int n_points = points.size();


    for (unsigned int i=0; i<n_points; ++i)
      {
        const Point<dim> &p = points[i];
        // Real:
        value_list[i](0) = 0.0;
        value_list[i](1) = 0.0;
        value_list[i](2) = 0.0;
        value_list[i](3) = 0.0;
        value_list[i](4) = 0.0;
        value_list[i](5) = 0.0;
      }
  }
  // END POLYNOMIALTEST
  template class polynomialTest<3>;
  
}