#include <boost/python/numpy.hpp>
#include <boost/scoped_array.hpp>
#include <iostream>
#include <math.h>

namespace p = boost::python;
namespace np = boost::python::numpy;

/*
void fill1(double * array, int rows, int cols, int row_stride, int col_stride) {
    double * row_iter = array;
    double n = 0.0;
    for (int i = 0; i < rows; ++i, row_iter += row_stride) {
        double * col_iter = row_iter;
        for (int j = 0; j < cols; ++j, col_iter += col_stride) {
            *col_iter = ++n;
        }
    }
}
*/

void print_ndarray2(double * array, int x_n, int y_n, int z_n, int x_str, int y_str, int z_str,int i, int j, int k);
void print_ndarray(double * array, int x_n, int y_n, int z_n, int x_str, int y_str, int z_str);

void print_ndarray(double * array, int x_n, int y_n, int z_n, int x_str, int y_str, int z_str) {
  double * x_iter = array;
  for (int i = 0; i < x_n; ++i, x_iter += x_str) {
    double * y_iter = x_iter;
    for (int j = 0; j < y_n; ++j, y_iter += y_str) {
      double * z_iter = y_iter;
      for (int k = 0; k < z_n; ++k, z_iter += z_str) {
	std::cout<<i<<" "<<j<<" "<<k<<" "<<*z_iter<<std::endl;
	print_ndarray2(array, x_n, y_n, z_n, x_str, y_str, z_str, i, j, k);
      }
    }
  }
}

void print_ndarray2(double * array, int x_n, int y_n, int z_n, int x_str, int y_str, int z_str,int i, int j, int k) {
  //double * x_iter = array;
  //for (int i = 0; i < x_n; ++i, x_iter += x_str) {
  //double * y_iter = x_iter;
  //for (int j = 0; j < y_n; ++j, y_iter += y_str) {
  //double * z_iter = y_iter;
  //for (int k = 0; k < z_n; ++k, z_iter += z_str) {
  std::cout<<i<<" "<<j<<" "<<k<<" "<<array[i*x_str + j*y_str + k*z_str]<<std::endl;
  //}
  //}
  //}
}

//void log_gaussian(const auto (&x), const auto (&mean), const auto (&sigma), auto (&log_pdf)){
void log_gaussian(np::ndarray const (&x),
		  np::ndarray const (&mean),
		  np::ndarray const (&sigma),
		  np::ndarray const (&log_pdf)){
  double *xd = reinterpret_cast<double*>(x.get_data());
  double *meand = reinterpret_cast<double*>(mean.get_data());
  double *sigmad = reinterpret_cast<double*>(sigma.get_data());
  double *log_pdfd = reinterpret_cast<double*>(log_pdf.get_data());
  const unsigned int xDim = mean.shape(0);
  const unsigned int yDim = mean.shape(1);
  const unsigned int zDim = mean.shape(2);
  const unsigned int x_str01 = x.strides(0) / sizeof(double);
  const unsigned int x_str02 = x.strides(1) / sizeof(double);
  const unsigned int x_str03 = x.strides(2) / sizeof(double);
  const unsigned int mean_str01 = mean.strides(0) / sizeof(double);
  const unsigned int mean_str02 = mean.strides(1) / sizeof(double);
  const unsigned int mean_str03 = mean.strides(2) / sizeof(double);
  const unsigned int sigma_str01 = sigma.strides(0) / sizeof(double);
  const unsigned int sigma_str02 = sigma.strides(1) / sizeof(double);
  const unsigned int sigma_str03 = sigma.strides(2) / sizeof(double);
  const unsigned int log_pdf_str01 = sigma.strides(0) / sizeof(double);
  const unsigned int log_pdf_str02 = sigma.strides(1) / sizeof(double);
  const unsigned int log_pdf_str03 = sigma.strides(2) / sizeof(double);
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  unsigned int i_x = 0;
  unsigned int i_mean = 0;
  unsigned int i_sigma = 0;
  unsigned int i_log_pdf = 0;
  //std::cout<<"xDim "<<xDim<<std::endl
  //	   <<"yDim "<<yDim<<std::endl
  //	   <<"zDim "<<zDim<<std::endl;
  for(i = 0;i<xDim;i++){
    for(j = 0;j<yDim;j++){
      for(k = 0;k<zDim;k++){
	i_x = i*x_str01 + j*x_str02 + k*x_str03;
	i_mean = i*mean_str01 + j*mean_str02 + k*mean_str03;
	i_sigma = i*sigma_str01 + j*sigma_str02 + k*sigma_str03;
	i_log_pdf = i*log_pdf_str01 + j*log_pdf_str02 + k*log_pdf_str03;
        log_pdfd[i_log_pdf] = -(xd[i_x] - meand[i_mean])*(xd[i_x] - meand[i_mean])/2.0/sigmad[i_sigma]/sigmad[i_sigma] - log(sqrt(2*3.141592653589793)*sigmad[i_sigma]);
	//log_pdfd[i_log_pdf] = -(xd[i_x] - meand[i_mean])*(xd[i_x] - meand[i_mean])/2.0/sigmad[i_sigma]/sigmad[i_sigma];
      }
    }
  }
}

void wrap_log_gaussian(np::ndarray const (&x),
		       np::ndarray const (&mean),
		       np::ndarray const (&sigma),
		       np::ndarray const (&log_pdf)) {
  log_gaussian(x,mean,sigma,log_pdf);
}

void wrap_print_ndarray(np::ndarray const & array) {
  /*
    if (array.get_dtype() != boost::python::numpy::dtype::get_builtin<double>()) {
        PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
        boost::python::throw_error_already_set();
    }
    if (array.get_nd() != 2) {
        PyErr_SetString(PyExc_TypeError, "Incorrect number of dimensions");
        boost::python::throw_error_already_set();
    }
  */
    std::cout<<"array.get_nd()  = "<<array.get_nd()<<std::endl;
    /*
    std::cout<<"sizeof(char)  = "<<sizeof(char)<<std::endl;
    std::cout<<"sizeof(int)  = "<<sizeof(int)<<std::endl;
    std::cout<<"sizeof(float)  = "<<sizeof(float)<<std::endl;
    std::cout<<"sizeof(double)  = "<<sizeof(double)<<std::endl;
    */
    
    for(int i = 0 ;i<array.get_nd();i++)
      std::cout<<"array.shape   = "<<array.shape(i)<<std::endl;
    for(int i = 0 ;i<array.get_nd();i++)
      std::cout<<"array.strides = "<<array.strides(i)/sizeof(double)<<std::endl;

    print_ndarray(reinterpret_cast<double*>(array.get_data()),
		  array.shape(0), array.shape(1), array.shape(2),
		  array.strides(0) / sizeof(double), array.strides(1) / sizeof(double), array.strides(2) / sizeof(double));
    
    
    /*
    for(int i = 0 ;i<array.get_nd();i++)
      std::cout<<"array.strides = "<<array.strides(i)<<std::endl;
    for(int i = 0 ;i<array.get_nd();i++)
    */

    //std::cout<<"array.shape(2) = "<<array.shape(2)<<std::endl;
    //fill1(reinterpret_cast<double*>(array.get_data()),
    //    array.shape(0), array.shape(1),
    //    array.strides(0) / sizeof(double), array.strides(1) / sizeof(double));
}

BOOST_PYTHON_MODULE(log_gaussian_boost_python) {
    np::initialize();
    p::def("print_ndarray", wrap_print_ndarray);
    p::def("log_gaussian", wrap_log_gaussian);
}
