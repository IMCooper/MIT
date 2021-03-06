#include <deal.II/base/polynomial.h>
#include <deal.II/base/point.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/thread_management.h>

#include <cmath>
#include <algorithm>
#include <limits>

#ifndef mypolynomials_h
#define mypolynomials_h

namespace MyPolynomials
{
  using namespace dealii;
  using namespace Polynomials;

  // Copy of the original deal.II legendre class, but with the coefficents adjusted
  // so that the recursive formula is for the integrated legendre polynomials described
  // in the thesis of Sabine Zaglmayer.
  //
  // L0 = -1  (added so that it can be generated recursively from 0)
  // L1 = x
  // L2 = 0.5*(x^2 - 1)
  // (n+1)Lnp1 = (2n-1)xLn - (n-2)Lnm1.
  //
  // However, it is possible to generate them directly from the legendre polynomials:
  //
  // (2n-1)Ln = ln - lnm2
  //
  // Could simply call the legendre set and adjust accordingly??
  class integratedLegendre : public Polynomial<double>
  {
  public:
    integratedLegendre (unsigned int p);
    
    static
    std::vector<Polynomial<double> >
    generate_complete_basis (const unsigned int degree); 
    
  private:
    static std::vector<std_cxx11::shared_ptr<const std::vector<double> > > shifted_coefficients;
    
    static std::vector<std_cxx11::shared_ptr<const std::vector<double> > > recursive_coefficients;

    static void compute_coefficients (const unsigned int p);
    
    static const std::vector<double> &
    get_coefficients (const unsigned int k);
    
  };
  namespace
  {
  Threads::Mutex coefficients_lock;
  }


// Reserve space for polynomials up to degree 19. Should be sufficient
// for the start.
  std::vector<std_cxx11::shared_ptr<const std::vector<double> > >
  integratedLegendre::recursive_coefficients(20);
  std::vector<std_cxx11::shared_ptr<const std::vector<double> > >
  integratedLegendre::shifted_coefficients(20);


  integratedLegendre::integratedLegendre (const unsigned int k)
    :
    Polynomial<double> (get_coefficients(k))
  {}

  void
  integratedLegendre::compute_coefficients (const unsigned int k_)
  {
    // NOT REQUIRED AT THE MOMENT:
    
    // make sure we call the
    // Polynomial::shift function
    // only with an argument with
    // which it will not crash the
    // compiler
// #ifdef DEAL_II_LONG_DOUBLE_LOOP_BUG
//     typedef double SHIFT_TYPE;
// #else
//     typedef long double SHIFT_TYPE;
// #endif

    unsigned int k = k_;

    // first make sure that no other
    // thread intercepts the
    // operation of this function;
    // for this, acquire the lock
    // until we quit this function
    Threads::Mutex::ScopedLock lock(coefficients_lock);

    // The first 2 coefficients are hard-coded
    if (k==0)
      k=1;
    // check: does the information
    // already exist?
    if ((recursive_coefficients.size() < k+1) ||
        ((recursive_coefficients.size() >= k+1) &&
         (recursive_coefficients[k] ==
          std_cxx11::shared_ptr<const std::vector<double> >())))
      // no, then generate the
      // respective coefficients
      {
        // make sure that there is enough space in the array for the coefficients,
        // so we have to resize it to size k+1

        // but it's more complicated than that: we call this function recursively, so if we simply
        // resize it to k+1 here, then compute the coefficients for degree k-1 by calling this
        // function recursively, then it will reset the size to k -- not enough for what we want to do below. the
        // solution therefore is to only resize the size if we are going to *increase* it
        if (recursive_coefficients.size() < k+1)
        {
          recursive_coefficients.resize (k+1);
        }
        if (k<=1)
        {
            // create coefficients vectors for k=0 and k=1
            //
            // allocate the respective later assign it to the coefficients array to make it const
            std::vector<double> *c0 = new std::vector<double>(1);
            (*c0)[0] = -1.;

            std::vector<double> *c1 = new std::vector<double>(2);
            (*c1)[0] = 0.;
            (*c1)[1] = 1.;

            // now make these arrays const. use shared_ptr for recursive_coefficients because
            // that avoids a memory leak that would appear if we used plain pointers.
            recursive_coefficients[0] =
              std_cxx11::shared_ptr<const std::vector<double> >(c0);
            recursive_coefficients[1] =
              std_cxx11::shared_ptr<const std::vector<double> >(c1);

            // Compute polynomials orthogonal on [0,1]
            /* REMOVED: NO SHIFT FOR NOW:
            c0 = new std::vector<double>(*c0);
            c1 = new std::vector<double>(*c1);
                       
            Polynomial<double>::shift<SHIFT_TYPE> (*c0, -1.);
            Polynomial<double>::scale(*c0, 2.);
            Polynomial<double>::shift<SHIFT_TYPE> (*c1, -1.);
            Polynomial<double>::scale(*c1, 2.);
            Polynomial<double>::multiply(*c1, std::sqrt(3.));

            shifted_coefficients[0]=std_cxx11::shared_ptr<const std::vector<double> >(c0);
            shifted_coefficients[1]=std_cxx11::shared_ptr<const std::vector<double> >(c1);
            */
          }
        else
          {
            // for larger numbers, compute the coefficients recursively. to do so, we have to release the
            // lock temporarily to allow the called function to acquire it itself
            coefficients_lock.release ();
            compute_coefficients(k-1);
            coefficients_lock.acquire ();

            std::vector<double> *ck = new std::vector<double>(k+1);

            
            const double a = 1.0/(k);
            const double b = a*(2.0*k-3.0);
            const double c = a*(k-3.0);

            (*ck)[k]   = b*(*recursive_coefficients[k-1])[k-1];
            (*ck)[k-1] = b*(*recursive_coefficients[k-1])[k-2];
            for (unsigned int i=1; i<= k-2 ; ++i)
            {
              (*ck)[i] = b*(*recursive_coefficients[k-1])[i-1]
                         -c*(*recursive_coefficients[k-2])[i];
            }
            (*ck)[0]   = -c*(*recursive_coefficients[k-2])[0];

            // finally assign the newly created vector to the const pointer in the/ coefficients array
            recursive_coefficients[k] =
              std_cxx11::shared_ptr<const std::vector<double> >(ck);
              
            // and compute the coefficients for [0,1]
            /* REMOVED THE SHIFT FOR NOW.
            ck = new std::vector<double>(*ck);
            
            Polynomial<double>::shift<SHIFT_TYPE> (*ck, -1.);
            Polynomial<double>::scale(*ck, 2.);
            Polynomial<double>::multiply(*ck, std::sqrt(2.*k+1.));
            
            shifted_coefficients[k] =
              std_cxx11::shared_ptr<const std::vector<double> >(ck);
            */
          };
      };
  }



  const std::vector<double> &
  integratedLegendre::get_coefficients (const unsigned int k)
  {
    // first make sure the coefficients get computed if so necessary
    compute_coefficients (k);

    // then get a pointer to the array of coefficients. do that in a MT safe way
    Threads::Mutex::ScopedLock lock (coefficients_lock);
    return *recursive_coefficients[k];
//     return *shifted_coefficients[k];
  }


  std::vector<Polynomial<double> >
  integratedLegendre::generate_complete_basis (const unsigned int degree)
  {
    std::vector<Polynomial<double> > v;
    v.reserve(degree+1);
    for (unsigned int i=0; i<=degree; ++i)
    {
      v.push_back (integratedLegendre(i));
    }
    return v;
  }
}

#endif