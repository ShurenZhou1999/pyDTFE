#include <vector>
#include <algorithm>
#include <math.h>
#include <iterator>
#include <boost/multi_array.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h> /* CGAL kernel */
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>  /* periodic DT */
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>         /* periodic DT */
#include <CGAL/Delaunay_triangulation_3.h>                    /* regular DT  */
#include <CGAL/natural_neighbor_coordinates_3.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Triangle_3.h>


// Traits and triangulation data structures
typedef CGAL::Exact_predicates_exact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_with_info_3<size_t,K>   VertexInfo;
typedef CGAL::Triangulation_data_structure_3<VertexInfo>        Tds3;
typedef CGAL::Delaunay_triangulation_3<K, Tds3>                 Delaunay;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>     P3Traits;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<P3Traits>     PDelaunay;


typedef Delaunay::Point                 Point;
typedef PDelaunay::Point                PPoint;
typedef K::FT                           FT;
typedef PDelaunay::Iso_cuboid           Iso_cuboid;

// Geometric object types
typedef CGAL::Point_3<K>                Point_3;
typedef CGAL::Sphere_3<K>               Sphere_3;
typedef CGAL::Tetrahedron_3<K>          Tetrahedron_3;
typedef CGAL::Triangle_3<K>             Triangle_3;

#define NDIM 3

using namespace std;



// ----------------------------------------------------------------------------------------

class DToutput
{
    public:
    vector<double> volume;    // Volume of the tetrahedron
    vector<double> vertex_volume;   // Volume summation of the tetrahedron around the vertex
    vector<double> field_rho;       // density esitmated at the vertex
    vector< vector<double> > field_pos;    // position at the vertex
    vector< vector<double> > field_vel;    // velocity at the vertex
    size_t Nsimplices;
    Delaunay tess;
    bool periodic;
};


template <typename T>
class MeshOutput_3D
{
    public:
    boost::multi_array<T, 3> Mesh_rho;
    boost::multi_array<T, 3> Mesh_vx;
    boost::multi_array<T, 3> Mesh_vy;
    boost::multi_array<T, 3> Mesh_vz;
    MeshOutput_3D(int Ngrid[NDIM])
    {
        Mesh_rho.resize(boost::extents[Ngrid[0]][Ngrid[1]][Ngrid[2]]); 
        Mesh_vx.resize(boost::extents[Ngrid[0]][Ngrid[1]][Ngrid[2]]);
        Mesh_vy.resize(boost::extents[Ngrid[0]][Ngrid[1]][Ngrid[2]]);
        Mesh_vz.resize(boost::extents[Ngrid[0]][Ngrid[1]][Ngrid[2]]);
    }
};



// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------

/* pad the 3D Vector field */
inline void insert_array3D( 
    vector<vector<double>> &A, 
    vector<vector<double>> &B, 
    size_t index, 
    double a0, double a1, double a2  )
{
    A.push_back( vector<double>{ A[index][0]+a0, A[index][1]+a1, A[index][2]+a2 } );
    B.push_back( B[index] );
}

/* pad the 3D scalar field */
inline void insert_array3D( 
    vector<vector<double>> &A, 
    vector<double> &B, 
    size_t index, 
    double a0, double a1, double a2  )
{
    A.push_back( vector<double>{ A[index][0]+a0, A[index][1]+a1, A[index][2]+a2 } );
    B.push_back( B[index] );
}


/* pad the velocity/scalar catalog */
template <typename VT>
void PaddingData( 
    vector<vector<double>>& XYZ , 
    VT& PaddedField , 
    double L[NDIM] = nullptr,
    double PaddingRate = 0.03
    )
{
    size_t Nsize = XYZ.size();

    double Lx = L[0], Ly = L[1], Lz = L[2];
    double ddx = PaddingRate * Lx, ddy = PaddingRate * Ly, ddz = PaddingRate * Lz;
    double Lx_min = 0 +ddx, Ly_min = 0 +ddy, Lz_min = 0 +ddz;
    double Lx_max = Lx-ddx, Ly_max = Ly-ddy, Lz_max = Lz-ddz;
    for(size_t i=0; i<Nsize; i++)
    {
        if( XYZ[i][0] < Lx_min ) insert_array3D( XYZ, PaddedField, i, Lx, 0, 0 );
        else if( XYZ[i][0] > Lx_max ) insert_array3D( XYZ, PaddedField, i, -Lx, 0, 0 );
        // ---------------------------------------------------------------------
        if( XYZ[i][1] < Ly_min )
        {
            insert_array3D( XYZ, PaddedField, i, 0, Ly, 0 );
            if( XYZ[i][0] < Lx_min ) insert_array3D( XYZ, PaddedField, i, Lx, Ly, 0 );
            else if( XYZ[i][0] > Lx_max ) insert_array3D( XYZ, PaddedField, i, -Lx, Ly, 0 );
        }
        // ---------------------------------------------------------------------
        else if( XYZ[i][1] > Ly_max )
        {
            insert_array3D( XYZ, PaddedField, i, 0, -Ly, 0 );
            if( XYZ[i][0] < Lx_min ) insert_array3D( XYZ, PaddedField, i, Lx, -Ly, 0 );
            else if( XYZ[i][0] > Lx_max ) insert_array3D( XYZ, PaddedField, i, -Lx, -Ly, 0 );
        }
        // ---------------------------------------------------------------------
        if( XYZ[i][2] < Lz_min )
        {
            insert_array3D( XYZ, PaddedField, i, 0, 0, Lz );
            if( XYZ[i][0] < Lx_min ) insert_array3D( XYZ, PaddedField, i, Lx, 0, Lz );
            else if( XYZ[i][0] > Lx_max ) insert_array3D( XYZ, PaddedField, i, -Lx, 0, Lz );
            if( XYZ[i][1] < Ly_min )
            {
                insert_array3D( XYZ, PaddedField, i, 0, Ly, Lz );
                if( XYZ[i][0] < Lx_min ) insert_array3D( XYZ, PaddedField, i, Lx, Ly, Lz );
                else if( XYZ[i][0] > Lx_max ) insert_array3D( XYZ, PaddedField, i, -Lx, Ly, Lz );
            }
            else if( XYZ[i][1] > Ly_max )
            {
                insert_array3D( XYZ, PaddedField, i, 0, -Ly, Lz );
                if( XYZ[i][0] < Lx_min ) insert_array3D( XYZ, PaddedField, i, Lx, -Ly, Lz );
                else if( XYZ[i][0] > Lx_max ) insert_array3D( XYZ, PaddedField, i, -Lx, -Ly, Lz );
            }
        }
        // ---------------------------------------------------------------------
        else if( XYZ[i][2] > Lz_max )
        {
            insert_array3D( XYZ, PaddedField, i, 0, 0, -Lz );
            if( XYZ[i][0] < Lx_min ) insert_array3D( XYZ, PaddedField, i, Lx, 0, -Lz );
            else if( XYZ[i][0] > Lx_max ) insert_array3D( XYZ, PaddedField, i, -Lx, 0, -Lz );
            if( XYZ[i][1] < Ly_min )
            {
                insert_array3D( XYZ, PaddedField, i, 0, Ly, -Lz );
                if( XYZ[i][0] < Lx_min ) insert_array3D( XYZ, PaddedField, i, Lx, Ly, -Lz );
                else if( XYZ[i][0] > Lx_max ) insert_array3D( XYZ, PaddedField, i, -Lx, Ly, -Lz );
            }
            else if( XYZ[i][1] > Ly_max )
            {
                insert_array3D( XYZ, PaddedField, i, 0, -Ly, -Lz );
                if( XYZ[i][0] < Lx_min ) insert_array3D( XYZ, PaddedField, i, Lx, -Ly, -Lz );
                else if( XYZ[i][0] > Lx_max ) insert_array3D( XYZ, PaddedField, i, -Lx, -Ly, -Lz );
            }
        }
        // ---------------------------------------------------------------------
    }
}




// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------


DToutput To_Delaunay_3D( 
            vector<vector<double>> XYZ 
            // double vector<double> Weight = nullptr
            )
{
    DToutput output;

    size_t Nsize = XYZ.size();
    vector< pair<Point,size_t> > points;
    for(size_t i=0; i<Nsize; i++)
        points.push_back( make_pair( Point_3(XYZ[i][0], XYZ[i][1], XYZ[i][2]), i ) );
    
    Delaunay tess(points.begin(), points.end());
    assert(tess.is_valid());
    if(points.size()!=tess.number_of_vertices())
        P_WRN <<"==> Warning: Number of vertices does not match the input particle numebr! (Possible repeated particles?)"<<std::endl;
    points.clear();
    /* tess.number_of_cells() --> includes the cell with infinite vertex, which is artificially assigned for dealing the boundary */
    /* tess.number_of_finite_cells() --> only counts the cells whose all 4 vertex are finite */
    size_t Nsimp = tess.number_of_finite_cells();
    size_t Nvert = tess.number_of_vertices();
    output.Nsimplices = Nsimp;
    
    output.vertex_volume.reserve(Nvert);
    for(size_t i=0; i<Nvert; output.vertex_volume[i++]=0 );
    
    Point_3 ipoint;
    Tetrahedron_3 itetra;
    int icell = 0;
    
    for(Delaunay::Finite_cells_iterator cell=tess.finite_cells_begin() ; cell!=tess.finite_cells_end() ; cell++ )
    {
        itetra = Tetrahedron_3( cell->vertex(0)->point(),
                                cell->vertex(1)->point(),
                                cell->vertex(2)->point(),
                                cell->vertex(3)->point() ) ;
        double ivolume = CGAL::to_double(itetra.volume());
        for(size_t i=0;i<4;i++)
            output.vertex_volume[cell->vertex(i)->info()] += ivolume;
        icell++;
    }
    if(icell-Nsimp!=0)
        P_WRN <<"==> Warning: Number of simplices does not match."<<std::endl;
    output.field_rho.resize(Nvert);

    // If loop over the vertex
    //for( Delaunay::Vertex_handle vertex = tess.vertices_begin(); vertex != tess.vertices_end(); vertex++ )
    // { size_t ivert = vertex->info();
    for(size_t ivert=0; ivert<Nvert; ivert++ )
        output.field_rho[ ivert ] = 4.0 / output.vertex_volume[ivert];   /* The density does not count the particle mass */
    output.tess = tess;

    return output;
}



// ----------------------------------------------------------------------------------------
/* Interpolate the density & velocity field on 3D grid mesh */
// ----------------------------------------------------------------------------------------

MeshOutput_3D<float> InterploateGrid_3D( 
            DToutput output,
            double L_min[NDIM], double L_max[NDIM], 
            int Ngrid[NDIM]
            )
{
    vector<double> dx = { (L_max[0]-L_min[0])/Ngrid[0], 
                          (L_max[1]-L_min[1])/Ngrid[1], 
                          (L_max[2]-L_min[2])/Ngrid[2] };
    MeshOutput_3D<float> meshValues(Ngrid);

    Delaunay dt = output.tess;
    Delaunay::Locate_type lt;
    Delaunay::Cell_handle current;

    int li, lj;
    double Gridx, Gridy, Gridz;
    Gridx = L_min[0] -0.5* dx[0];
    for( int i0=0 ; i0<Ngrid[0]; i0++ )
    {
        Gridx += dx[0];
        Gridy = L_min[1] -0.5* dx[1];
        for( int i1=0 ; i1<Ngrid[1]; i1++ )
        {
            Gridy += dx[1];
            Gridz = L_min[2] -0.5* dx[2];
            for( int i2=0 ; i2<Ngrid[2]; i2++ )
            {
                Gridz += dx[2];
                Point isample(Gridx, Gridy, Gridz);
                current = dt.locate(isample, lt, li, lj) ;
                /* ------------------------------------------------------------------------------ */

                double diff_den[NDIM];             // vector, $\rho_{\li} - \rho_{l0}$
                double diff_pos[NDIM][NDIM];       // matrix, in i-th row $\vec{x}_{li} - \vec{x}_{l0}$
                double diff_vel[NDIM][NDIM];       // matrix, in i-th row $\vec{v}_{li} - \vec{v}_{l0}$
                size_t ivert0 = current->vertex(0)->info() ;    // base vertex
                Point baseP  = current->vertex(0)->point();
                for(int i=0; i<NDIM; i++)     // over 3 vertex, other than the base
                {
                    size_t ivert = current->vertex(i+1)->info() ;
                    diff_den[i] = output.field_rho[ivert] - output.field_rho[ivert0] ;
                    for(int j=0; j<NDIM; j++)    // over 3 dimensions (x, y, z)
                    {
                        diff_pos[i][j] = CGAL::to_double( current->vertex(i+1)->point()[j] - baseP[j] );
                        diff_vel[i][j] = output.field_vel[ivert][j] - output.field_vel[ivert0][j] ;    // The velocity is stored as (dims, particles)
                    }
                }
                double posMatrixInverse[NDIM][NDIM];
                double Grad_den[NDIM];
                double Grad_vel[NDIM][NDIM];
                matrixInverse_3x3( diff_pos, posMatrixInverse );
                matrixMultiplication( posMatrixInverse, diff_den, Grad_den );    //computes the density gradient = posMatrixInverse * dens
                matrixMultiplication( posMatrixInverse, diff_vel, Grad_vel );
                matrixTranspose(Grad_vel);
                
                double dpos[NDIM] = {Gridx - CGAL::to_double(baseP[0]), 
                                     Gridy - CGAL::to_double(baseP[1]), 
                                     Gridz - CGAL::to_double(baseP[2]) };
                meshValues.Mesh_rho[i0][i1][i2] = (float)( output.field_rho[ivert0]    + dotProduct3(Grad_den, dpos) );
                meshValues.Mesh_vx[i0][i1][i2]  = (float)( output.field_vel[ivert0][0] + dotProduct3(Grad_vel[0], dpos) );
                meshValues.Mesh_vy[i0][i1][i2]  = (float)( output.field_vel[ivert0][1] + dotProduct3(Grad_vel[1], dpos) );
                meshValues.Mesh_vz[i0][i1][i2]  = (float)( output.field_vel[ivert0][2] + dotProduct3(Grad_vel[2], dpos) );
                /* ------------------------------------------------------------------------------ */
            }
        }
    }
    return meshValues;
}



// ----------------------------------------------------------------------------------------
/* Interpolate the density & velocity field on given sampling points */
// ----------------------------------------------------------------------------------------

vector<vector<float>> Interploate_3D( 
            DToutput output,
            vector<vector<double>> samplingPoint
            )
{
    Delaunay dt = output.tess;
    Delaunay::Locate_type lt;
    int li, lj;
    Delaunay::Cell_handle current;

    int Nsamples = samplingPoint[0].size();
    vector<vector<float>> meshValues(4, vector<float>(Nsamples, 0.));  // (Nkind, Nsamples)

    for(int isam=0; isam<Nsamples; isam++)
    {
        double Gridx = samplingPoint[0][isam];
        double Gridy = samplingPoint[1][isam];
        double Gridz = samplingPoint[2][isam];
        Point isample(Gridx, Gridy, Gridz);
        current = dt.locate(isample, lt, li, lj) ;
        /* ------------------------------------------------------------------------------ */

        double diff_den[NDIM];             // vector, $\rho_{\li} - \rho_{l0}$
        double diff_pos[NDIM][NDIM];       // matrix, in i-th row $\vec{x}_{li} - \vec{x}_{l0}$
        double diff_vel[NDIM][NDIM];       // matrix, in i-th row $\vec{v}_{li} - \vec{v}_{l0}$
        size_t ivert0 = current->vertex(0)->info() ;    // base vertex
        Point baseP  = current->vertex(0)->point();
        for(int i=0; i<NDIM; i++)     // over 3 vertex, other than the base
        {
            size_t ivert = current->vertex(i+1)->info() ;
            diff_den[i] = output.field_rho[ivert] - output.field_rho[ivert0] ;
            for(int j=0; j<NDIM; j++)    // over 3 dimensions (x, y, z)
            {
                diff_pos[i][j] = CGAL::to_double( current->vertex(i+1)->point()[j] - baseP[j] );
                diff_vel[i][j] = output.field_vel[ivert][j] - output.field_vel[ivert0][j] ;    // The velocity is stored as (dims, particles)
            }
        }
        double posMatrixInverse[NDIM][NDIM];
        double Grad_den[NDIM];
        double Grad_vel[NDIM][NDIM];
        matrixInverse_3x3( diff_pos, posMatrixInverse );
        matrixMultiplication( posMatrixInverse, diff_den, Grad_den );    //computes the density gradient = posMatrixInverse * dens
        matrixMultiplication( posMatrixInverse, diff_vel, Grad_vel );
        matrixTranspose(Grad_vel);
        
        double dpos[NDIM] ={Gridx - CGAL::to_double(baseP[0]), 
                            Gridy - CGAL::to_double(baseP[1]), 
                            Gridz - CGAL::to_double(baseP[2]) };
        meshValues[0][isam] = (float)( output.field_rho[ivert0]    + dotProduct3(Grad_den, dpos) );
        meshValues[1][isam] = (float)( output.field_vel[ivert0][0] + dotProduct3(Grad_vel[0], dpos) );
        meshValues[2][isam] = (float)( output.field_vel[ivert0][1] + dotProduct3(Grad_vel[1], dpos) );
        meshValues[3][isam] = (float)( output.field_vel[ivert0][2] + dotProduct3(Grad_vel[2], dpos) );
        /* ------------------------------------------------------------------------------ */
    }
    return meshValues;
}



// ----------------------------------------------------------------------------------------
/* Interpolate the scalar field on given sampling points */
// ----------------------------------------------------------------------------------------

vector<float> Interploate_3D_ScalarField( 
            DToutput output,
            vector<vector<double>> samplingPoint, 
            vector<double> scalarField       // same size as particle position number
            )
{
    Delaunay dt = output.tess;
    Delaunay::Locate_type lt;
    int li, lj;
    Delaunay::Cell_handle current;

    int Nsamples = samplingPoint[0].size();
    vector<float> out_dens;
    vector<float> out_scal;

    for(int isam=0; isam<Nsamples; isam++)
    {
        double Gridx = samplingPoint[0][isam];
        double Gridy = samplingPoint[1][isam];
        double Gridz = samplingPoint[2][isam];
        Point isample(Gridx, Gridy, Gridz);
        current = dt.locate(isample, lt, li, lj) ;
        /* ------------------------------------------------------------------------------ */

        double diff_den[NDIM];
        double diff_pos[NDIM][NDIM];
        double diff_scal[NDIM];
        size_t ivert0 = current->vertex(0)->info() ;
        Point baseP  = current->vertex(0)->point();
        for(int i=0; i<NDIM; i++)
        {
            size_t ivert = current->vertex(i+1)->info() ;
            diff_den[i]  = output.field_rho[ivert] - output.field_rho[ivert0] ;
            diff_scal[i] = scalarField[ivert] - scalarField[ivert0] ;
            for(int j=0; j<NDIM; j++)
                diff_pos[i][j] = CGAL::to_double( current->vertex(i+1)->point()[j] - baseP[j] );
        }
        double posMatrixInverse[NDIM][NDIM];
        double Grad_den[NDIM];
        double Grad_scal[NDIM];
        matrixInverse_3x3( diff_pos, posMatrixInverse );
        matrixMultiplication( posMatrixInverse, diff_den,  Grad_den  );    //computes the density gradient = posMatrixInverse * dens
        matrixMultiplication( posMatrixInverse, diff_scal, Grad_scal );
        
        double dpos[NDIM] ={Gridx - CGAL::to_double(baseP[0]), 
                            Gridy - CGAL::to_double(baseP[1]), 
                            Gridz - CGAL::to_double(baseP[2]) };
        out_dens[isam] = (float)( output.field_rho[ivert0] + dotProduct3(Grad_den , dpos) );
        out_scal[isam] = (float)( scalarField[ivert0]      + dotProduct3(Grad_scal, dpos) );
        /* ------------------------------------------------------------------------------ */
    }
    return out_scal;
}

