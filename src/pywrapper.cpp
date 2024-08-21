#include<iostream>
#include<math.h>
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>

namespace py = pybind11;
using namespace std;
using namespace py::literals;



// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------


vector<vector<double>> NdataToVector( py::array_t<float> arr_in_ )
{
    /*
        Convert Numpy Array to C++ data type: vector<vector<double>>
    */
    int Ndim = arr_in_.shape(1) ;      // Number of dimensions (2D/3D particle data)
    auto Np  = arr_in_.shape(0) ;
    if(Ndim != 2 && Ndim != 3)
        throw runtime_error("NdataToVector::Input data must be 2D or 3D.");
    auto arr_in = arr_in_.mutable_unchecked<2>();
    vector<vector<double>> arr_out( Ndim, vector<double>(Np, 0.0) );
    for(int i=0; i<Ndim; i++)
        for(int j=0; j<Np; j++)
            arr_out[i][j] = static_cast<double>(arr_in(j, i));
    return arr_out;
}


vector<double> VectorToVector( py::array_t<float> arr_in_ )
{
    py::buffer_info buffer = arr_in_.request();
    if( buffer.ndim != 1 )
        throw runtime_error("VectorToVector::Input data must be 1D.");
    auto Nd     = arr_in_.shape(0) ;
    auto arr_in = arr_in_.mutable_unchecked<1>();
    vector<double> arr_out(Nd, 0.0);
    for(int i=0; i<Nd; arr_out[++i]=static_cast<double>(arr_in(i)) );
    return arr_out;
}




// ----------------------------------------------------------------------------------------
/* Interpolation on uniform grid */
// ----------------------------------------------------------------------------------------

py::array_t<double> DTFE_3D_Grid( 
            py::array_t<float> arr_pos_, 
            py::array_t<float> arr_vel_, 
            int Nmesh , 
            double Boxsize, 
            double paddingRate=0.05
            )
{
    // py::print( "DTFE_3D::running ... " );               ////////////////////////////////
    py::module_ np = py::module_::import("numpy");
    py::tuple args = py::make_tuple( 4, Nmesh, Nmesh, Nmesh );
    py::array_t<double> arr_out_ = np.attr("zeros")( "shape"_a = args );
    auto arr_out = arr_out_.mutable_unchecked<4>();
    auto Np = arr_pos_.shape(0) ;
    int dim = 3;

    vector<vector<double>> pos = NdataToVector( arr_pos_ );
    vector<vector<double>> vel = NdataToVector( arr_vel_ );

    double L[dim] = { Boxsize, Boxsize, Boxsize };
    int Ngrid[dim] = { Nmesh, Nmesh, Nmesh };
    // py::print( "==> Building Delaunay Triangulation." );                   ////////////////////////////////
    if(paddingRate>1e-10)
        PaddingData(pos, vel, L, paddingRate);
    DToutput output = To_Delaunay_3D( pos );
    output.field_vel = vel ;
    // py::print( "==> Constructing grid field ... " );                   ////////////////////////////////
    MeshOutput_3D<double> meshValues = InterploateGrid_3D( output, L, Ngrid );

    for(int i0=0; i0<Nmesh; i0++)
        for(int i1=0; i1<Nmesh; i1++)
            for(int i2=0; i2<Nmesh; i2++)
            {
                arr_out(0, i0, i1, i2) = meshValues.Mesh_rho[i0][i1][i2];
                arr_out(1, i0, i1, i2) = meshValues.Mesh_vx[ i0][i1][i2];
                arr_out(2, i0, i1, i2) = meshValues.Mesh_vy[ i0][i1][i2];
                arr_out(3, i0, i1, i2) = meshValues.Mesh_vz[ i0][i1][i2];
            }
    // py::print( "DTFE_3D::End." );               ////////////////////////////////
    return arr_out_;
}



py::array_t<double> DTFE_2D_Grid( 
            py::array_t<float> arr_pos_, 
            py::array_t<float> arr_vel_, 
            int Nmesh , 
            double Boxsize, 
            double paddingRate=0.05
            )
{
    // py::print( "DTFE_2D::running ... " );               ////////////////////////////////
    py::module_ np = py::module_::import("numpy");
    py::tuple args = py::make_tuple( 3, Nmesh, Nmesh );
    py::array_t<double> arr_out_ = np.attr("zeros")( "shape"_a = args );
    auto arr_out = arr_out_.mutable_unchecked<3>();

    vector<vector<double>> pos = NdataToVector( arr_pos_ );
    vector<vector<double>> vel = NdataToVector( arr_vel_ );
    auto Np = arr_pos_.shape(0) ;
    int dim = 2 ;

    double L[dim] = { Boxsize, Boxsize };
    int Ngrid[dim] = { Nmesh, Nmesh };
    // py::print( "==> Building Delaunay Triangulation." );                   ////////////////////////////////
    if(paddingRate>1e-10)
        PaddingData_2D(pos, vel, L, paddingRate);
    DToutput2 output = To_Delaunay_2D( pos );
    output.field_vel = vel ;
    // py::print( "==> Constructing grid field ... " );                   ////////////////////////////////
    MeshOutput_2D<double> meshValues = InterploateGrid_2D( output, L, Ngrid );

    for(int i0=0; i0<Nmesh; i0++)
        for(int i1=0; i1<Nmesh; i1++)
        {
            arr_out(0, i0, i1) = meshValues.Mesh_rho[i0][i1];
            arr_out(1, i0, i1) = meshValues.Mesh_vx[ i0][i1];
            arr_out(2, i0, i1) = meshValues.Mesh_vy[ i0][i1];
        }
    // py::print( "DTFE_2D::End." );               ////////////////////////////////
    return arr_out_;
}




// ----------------------------------------------------------------------------------------
/* Interpolation on sampling points */
// ----------------------------------------------------------------------------------------

/* For both 3D and 2D case */
py::array_t<double> DTFE_SampleVel( 
            py::array_t<float> arr_pos_, 
            py::array_t<float> arr_vel_, 
            py::array_t<float> arr_sampling_,
            double Boxsize, 
            double paddingRate=0.05
            )
{
    py::module_ np = py::module_::import("numpy");
    auto Nsamp  = arr_sampling_.shape(0) ;
    auto Nshape = arr_sampling_.shape(1) ;
    py::tuple args = py::make_tuple( Nshape+1, Nsamp );
    py::array_t<double> arr_out_ = np.attr("zeros")( "shape"_a = args );
    auto arr_out = arr_out_.mutable_unchecked<2>();

    vector<vector<double>> pos   = NdataToVector( arr_pos_ );
    vector<vector<double>> vel   = NdataToVector( arr_vel_ );
    vector<vector<double>> samps = NdataToVector( arr_sampling_ );

    vector<vector<double>> meshValues;
    if(Nshape==3)
    {
        double L[3] = { Boxsize, Boxsize, Boxsize };
        if(paddingRate>1e-10)
            PaddingData(pos, vel, L, paddingRate);
        DToutput output = To_Delaunay_3D( pos );
        output.field_vel = vel ;
        meshValues = Interploate_3D( output, samps );
    }
    else if(Nshape==2)
    {
        double L[2] = { Boxsize, Boxsize };
        if(paddingRate>1e-10)
            PaddingData_2D(pos, vel, L, paddingRate);
        DToutput2 output = To_Delaunay_2D( pos );
        output.field_vel = vel ;
        meshValues = Interploate_2D( output, samps );
    }
    else
        throw runtime_error("DTFE_SampleVel::Invalid shape of sampling points.");

    for(int i0=0; i0<Nshape; i0++)
        for(int i1=0; i1<Nsamp; i1++)
            arr_out(i0, i1) = meshValues[i0][i1];
    return arr_out_;
}


/* For both 3D and 2D case */
py::array_t<double> DTFE_SampleScalar( 
            py::array_t<float> arr_pos_, 
            py::array_t<float> arr_scalar_, 
            py::array_t<float> arr_sampling_,
            double Boxsize, 
            double paddingRate=0.05
            )
{
    py::module_ np = py::module_::import("numpy");
    auto Nsamp  = arr_sampling_.shape(0) ;
    auto Nshape = arr_sampling_.shape(1) ;
    py::tuple args = py::make_tuple( Nsamp );
    py::array_t<double> arr_out_ = np.attr("zeros")( "shape"_a = args );
    auto arr_out = arr_out_.mutable_unchecked<1>();

    vector<vector<double>> pos      = NdataToVector(  arr_pos_      );
    vector<double>         scalar   = VectorToVector( arr_scalar_   );
    vector<vector<double>> sampling = NdataToVector(  arr_sampling_ );

    vector<double> meshValues;
    if(Nshape==3)
    {
        double L[3] = { Boxsize, Boxsize, Boxsize };
        if(paddingRate>1e-10)
            PaddingData(pos, scalar, L, paddingRate);
        DToutput output = To_Delaunay_3D( pos );
        meshValues = Interploate_3D_ScalarField( output, sampling, scalar );
    }
    else if(Nshape==2)
    {
        double L[2] = { Boxsize, Boxsize };
        if(paddingRate>1e-10)
            PaddingData_2D(pos, scalar, L, paddingRate);
        DToutput2 output = To_Delaunay_2D( pos );
        meshValues = Interploate_2D_ScalarField( output, sampling, scalar );
    }
    else
        throw runtime_error("DTFE_SampleScalar::Invalid shape of sampling points.");

    for(int i=0; i<Nsamp; arr_out(i++)=meshValues[i] );
    return arr_out_;
}


