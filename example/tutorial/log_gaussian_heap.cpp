#include <math.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

char const* greet(){
   return "hello, world";
}

void log_gaussian(const auto (&x), const auto (&mean), const auto (&sigma), auto (&log_pdf)){
  const unsigned int xDim = 196;
  const unsigned int yDim = 609;
  const unsigned int zDim = 38;
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  for(i = 0;i<xDim;i++)
    for(j = 0;j<yDim;j++)
      for(k = 0;k<zDim;k++)
	log_pdf[i][j][k] = -(x[i][j][0] - mean[i][j][k])*(x[i][j][0] - mean[i][j][k])/2.0/sigma[i][j][k]/sigma[i][j][k] - log(sqrt(2*3.141592653589793)*sigma[i][j][k]);
}

BOOST_PYTHON_MODULE(log_gaussian_heap_ext){
    using namespace boost::python;
    def("greet", greet);
    def("log_gaussian", log_gaussian);
}
