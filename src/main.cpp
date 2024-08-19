#include <vector>
#include <algorithm>
#include <math.h>
#include <iterator>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "./DTFE.h"
#include "./MathFuncs.cpp"
#include "./FileIO.cpp"
// #include "./DataType.cpp"
#include "./DTFE_3D.cpp"
#include "./DTFE_2D.cpp"
#include "./pywrapper.cpp"

namespace py = pybind11;
using namespace std;
using namespace py::literals;


/* ----------------------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------------------- */


PYBIND11_MODULE( src_dtfe, m )
{
    m.doc() = "DTFE:pybind11 module";
    m.def("DTFE_3D_Grid", &DTFE_3D_Grid, 
            "3D DTFE Velocity field interpolation on regular grid", 
            "position"_a, "velocity"_a, "Nmesh"_a, "Boxsize"_a, "paddingRate"_a  );
    m.def("DTFE_2D_Grid", &DTFE_2D_Grid, 
            "2D DTFE Velocity field interpolation on regular grid", 
            "position"_a, "velocity"_a, "Nmesh"_a, "Boxsize"_a, "paddingRate"_a  );

    m.def("DTFE_3D_SampleVel", &DTFE_3D_SampleVel, 
            "3D DTFE Velocity field interpolation on given samping points", 
            "position"_a, "velocity"_a, "sampling"_a, "Boxsize"_a, "paddingRate"_a  );
    m.def("DTFE_3D_SampleScalar", &DTFE_3D_SampleScalar, 
            "3D DTFE Scalar filed interpolation on given samping points", 
            "position"_a, "scalar"_a, "sampling"_a, "Boxsize"_a, "paddingRate"_a  );

    //m.def("Test", &Test, "Test", 
    //        "position"_a, "velocity"_a, "Nmesh"_a, "Boxsize"_a  );
};

