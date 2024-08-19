#include <vector>
#include <algorithm>
#include <math.h>
#include <iterator>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h> /* CGAL kernel */
#include <CGAL/Delaunay_triangulation_2.h>                    /* regular DT  */
#include <CGAL/natural_neighbor_coordinates_2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangle_2.h>


// Traits and triangulation data structures
typedef CGAL::Exact_predicates_exact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_with_info_2<size_t,K>   VertexInfo2;
typedef CGAL::Triangulation_data_structure_2<VertexInfo2>       Tds2;
typedef CGAL::Delaunay_triangulation_2<K, Tds2>                 Delaunay2;


typedef Delaunay2::Point                Point2;
typedef K::FT                           FT;

// Geometric object types
typedef CGAL::Point_2<K>                Point_2;
typedef CGAL::Triangle_2<K>             Triangle_2;

#define NDIM2 2

using namespace std;

// #include "./DTFE.h"
// #include "./MathFuncs.cpp"
// #include "./FileIO.cpp"



// ----------------------------------------------------------------------------------------

class DToutput2
{
    public:
    vector<double> area;    // area of the tetrahedron
    vector<double> vertex_area;   // area summation of the tetrahedron around the vertex
    vector<double> field_rho;       // density esitmated at the vertex
    vector< vector<double> > field_vel;    // velocity at the vertex
    size_t Nsimplices;
    Delaunay2 tess;
    bool periodic;
};


template <typename T>
class MeshOutput_2D
{
    public:
    vector<vector<T>> Mesh_rho ; 
    vector<vector<T>> Mesh_vx ; 
    vector<vector<T>> Mesh_vy ; 
    MeshOutput_2D(int Ngrid[NDIM2])
    {
        Mesh_rho.resize(Ngrid[0], vector<T>(Ngrid[1], 0.)); 
        Mesh_vx.resize( Ngrid[0], vector<T>(Ngrid[1], 0.)); 
        Mesh_vy.resize( Ngrid[0], vector<T>(Ngrid[1], 0.)); 
    }
};



// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------


/* pad the 2D Vector field */
inline void insert_array2D( 
    vector<vector<double>> &A, 
    vector<vector<double>> &B, 
    size_t index, 
    double a0, double a1 )
{
    A[0].push_back( A[0][index] +a0 );
    A[1].push_back( A[1][index] +a1 );
    B[0].push_back( B[0][index] );
    B[1].push_back( B[1][index] );
}

/* pad the 2D scalar field */
inline void insert_array2D( 
    vector<vector<double>> &A, 
    vector<double> &B, 
    size_t index, 
    double a0, double a1 )
{
    A[0].push_back( A[0][index] +a0 );
    A[1].push_back( A[1][index] +a1 );
    B.push_back( B[index] );
}


/* pad the velocity/scalar catalog */
template <typename VT>
void PaddingData_2D( 
    vector<vector<double>>& XY , 
    VT& PaddedField , 
    double L[NDIM] = nullptr,
    double PaddingRate = 0.03
    )
{
    size_t Nsize = XY[0].size();

    double Lx = L[0], Ly = L[1];
    double ddx = PaddingRate * Lx, ddy = PaddingRate * Ly ;
    double Lx_min = 0 +ddx, Ly_min = 0 +ddy ;
    double Lx_max = Lx-ddx, Ly_max = Ly-ddy ;
    for(size_t i=0; i<Nsize; i++)
    {
        if( XY[0][i] < Lx_min ) insert_array2D( XY, PaddedField, i, Lx, 0 );
        else if( XY[0][i] > Lx_max ) insert_array2D( XY, PaddedField, i, -Lx, 0 );
        if( XY[1][i] < Ly_min )
        {
            insert_array2D( XY, PaddedField, i, 0, Ly );
            if( XY[0][i] < Lx_min ) insert_array2D( XY, PaddedField, i, Lx, Ly );
            else if( XY[0][i] > Lx_max ) insert_array2D( XY, PaddedField, i, -Lx, Ly );
        }
        else if( XY[1][i] > Ly_max )
        {
            insert_array2D( XY, PaddedField, i, 0, -Ly );
            if( XY[0][i] < Lx_min ) insert_array2D( XY, PaddedField, i, Lx, -Ly );
            else if( XY[0][i] > Lx_max ) insert_array2D( XY, PaddedField, i, -Lx, -Ly );
        }
    }
}



// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------


DToutput2 To_Delaunay_2D( 
            vector<vector<double>> XY , 
            bool Periodic = true, 
            double L[NDIM2] = nullptr,
            double PaddingRate = 0.01
            )
{
    DToutput2 output;
    output.periodic = Periodic;

    size_t Nsize = XY[0].size();
    vector< pair<Point2,size_t> > points;
    for(size_t i=0; i<Nsize; i++)
        points.push_back( make_pair( Point_2(XY[0][i], XY[1][i]),i ) );
    
    Delaunay2 tess(points.begin(), points.end());
    assert(tess.is_valid());
    if(points.size()!=tess.number_of_vertices())
        P_WRN <<"==> Warning: Number of vertices does not match the input particle numebr! (Possible repeated particles?)"<<std::endl;
    points.clear();
    size_t Nsimp = tess.number_of_faces();
    size_t Nvert = tess.number_of_vertices();
    output.Nsimplices = Nsimp;
    
    output.vertex_area.reserve(Nvert);
    for(size_t i=0; i<Nvert; output.vertex_area[i++]=0 );
    
    Point_2 ipoint;
    Triangle_2 itetra;
    int iface = 0;
    
    for(Delaunay2::Finite_faces_iterator face=tess.finite_faces_begin() ; face!=tess.finite_faces_end() ; face++ )
    {
        itetra = Triangle_2( face->vertex(0)->point(),
                            face->vertex(1)->point(),
                            face->vertex(2)->point() ) ;
        double iarea = CGAL::to_double(itetra.area());
        for(size_t i=0;i<3;i++)
            output.vertex_area[face->vertex(i)->info()] += iarea;
        iface++;
    }
    if(iface-Nsimp!=0)
        P_WRN <<"==> Warning: Number of simplices does not match."<<std::endl;
    output.field_rho.resize(Nvert);

    for(size_t ivert=0; ivert<Nvert; ivert++ )
        output.field_rho[ ivert ] = 3.0 / output.vertex_area[ivert];   /* The density does not count the particle mass */
    output.tess = tess;

    return output;
}



// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------


MeshOutput_2D<double> InterploateGrid_2D( 
            DToutput2 output,
            double L[NDIM2], int Ngrid[NDIM2] 
            )
{
    vector<double> dx = { L[0]/Ngrid[0], L[1]/Ngrid[1] };
    MeshOutput_2D<double> meshValues(Ngrid);

    Delaunay2 dt = output.tess;
    Delaunay2::Locate_type lt;
    Delaunay2::Face_handle current;

    int li, lj;
    double Gridx, Gridy;
    Gridx = -0.5* dx[0];
    for( int i0=0 ; i0<Ngrid[0]; i0++ )
    {
        Gridx += dx[0];
        Gridy = -0.5* dx[1];
        for( int i1=0 ; i1<Ngrid[1]; i1++ )
        {
            Gridy += dx[1];
            Point2 isample(Gridx, Gridy);
            current = dt.locate(isample) ;
            /* ------------------------------------------------------------------------------ */

            double diff_den[NDIM2];             // vector, $\rho_{\li} - \rho_{l0}$
            double diff_pos[NDIM2][NDIM2];       // matrix, in i-th row $\vec{x}_{li} - \vec{x}_{l0}$
            double diff_vel[NDIM2][NDIM2];       // matrix, in i-th row $\vec{v}_{li} - \vec{v}_{l0}$
            int ivert0 = current->vertex(0)->info() ;    // base vertex
            Point2 baseP  = current->vertex(0)->point();
            for(int i=0; i<NDIM2; i++)     // over 3 vertex, other than the base
            {
                int ivert = current->vertex(i+1)->info() ;
                diff_den[i] = output.field_rho[ivert] - output.field_rho[ivert0] ;
                for(int j=0; j<NDIM2; j++)    // over 3 dimensions (x, y, z)
                {
                    diff_pos[i][j] = CGAL::to_double( current->vertex(i+1)->point()[j] - baseP[j] );
                    diff_vel[i][j] = output.field_vel[j][ivert] - output.field_vel[j][ivert0] ;    // The velocity is stored as (dims, particles)
                }
            }
            double posMatrixInverse[NDIM2][NDIM2];
            double Grad_den[NDIM2];
            double Grad_vel[NDIM2][NDIM2];
            matrixInverse_2x2( diff_pos, posMatrixInverse );
            matrixMultiplication( posMatrixInverse, diff_den, Grad_den );    //computes the density gradient = posMatrixInverse * dens
            matrixMultiplication( posMatrixInverse, diff_vel, Grad_vel );
            matrixTranspose(Grad_vel);
            
            double dpos[NDIM2] = {  Gridx - CGAL::to_double(baseP[0]), 
                                    Gridy - CGAL::to_double(baseP[1]) };
            meshValues.Mesh_rho[i0][i1] = output.field_rho[ivert0]    + dotProduct2(Grad_den, dpos);
            meshValues.Mesh_vx[i0][i1]  = output.field_vel[0][ivert0] + dotProduct2(Grad_vel[0], dpos);
            meshValues.Mesh_vy[i0][i1]  = output.field_vel[1][ivert0] + dotProduct2(Grad_vel[1], dpos);
            /* ------------------------------------------------------------------------------ */
        }
    }
    return meshValues;
}



// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------

