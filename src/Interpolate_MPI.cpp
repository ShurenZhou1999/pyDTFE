#include<iostream>
#include<math.h>
#include<mpi.h>
#include<boost/multi_array.hpp>

using namespace std;


template <typename VT>
void SubDivideBox( 
    vector<vector<double>>& XYZ , 
    VT& PaddedField , 
    double subL_min[3], double subL_max[3], 
    double PaddingRate = 0.05
    )
{
    vector<vector<double>> subXYZ;
    VT subField;
    double edgex_min = subL_min[0] - PaddingRate*(subL_max[0]-subL_min[0]);
    double edgex_max = subL_max[0] + PaddingRate*(subL_max[0]-subL_min[0]);
    double edgey_min = subL_min[1] - PaddingRate*(subL_max[1]-subL_min[1]);
    double edgey_max = subL_max[1] + PaddingRate*(subL_max[1]-subL_min[1]);
    double edgez_min = subL_min[2] - PaddingRate*(subL_max[2]-subL_min[2]);
    double edgez_max = subL_max[2] + PaddingRate*(subL_max[2]-subL_min[2]);
    vector<vector<double>> XYZ_sub;
    VT PaddedField_sub;
    for(int i=0; i<XYZ.size(); i++)
    {
        if( XYZ[i][0] > edgex_min )
        if( XYZ[i][0] < edgex_max )
        if( XYZ[i][1] > edgey_min )
        if( XYZ[i][1] < edgey_max )
        if( XYZ[i][2] > edgez_min )
        if( XYZ[i][2] < edgez_max )
        {
            subXYZ.push_back( XYZ[i] );
            subField.push_back( PaddedField[i] );
        }
    }
    XYZ.clear();
    XYZ.shrink_to_fit();
    PaddedField.clear();
    PaddedField.shrink_to_fit();
    // Deep Copy the sub-data to full-data, to save memory.
    XYZ = subXYZ;
    PaddedField = subField;
}



boost::multi_array<float, 4> DTFE_3D_Grid_MPI( 
            vector<vector<double>> pos,
            vector<vector<double>> vel,
            int Nmesh , 
            double Boxsize, 
            double paddingRate=0.05
            )
{
    /********* Subdivde the processes *********/
    // MPI_Init(NULL, NULL); 
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int nsplitx = cbrt(world_size);
    int nsplitComm = nsplitx ;
    int nsplity = nsplitx, nsplitz = nsplitx;
    if( (nsplitx+1)*(nsplitx+1)*(nsplitx+1) <= world_size)
        { nsplitx += 1; nsplity+=1; nsplitz+=1; nsplitComm=nsplitx+1; }
    else if ( nsplitx*(nsplitx+1)*(nsplitx+1) <= world_size )
        { nsplitx += 1; nsplity+=1; nsplitComm=nsplitx*(nsplitx+1);  }
    else if ( nsplitx*nsplitx*(nsplitx+1) <= world_size )
        { nsplitx += 1; nsplitComm=nsplitx*(nsplitx+1); }

    int Nprocess = nsplitx*nsplity*nsplitz ;
    if(world_rank >= Nprocess)
    {
        for(int i=0; i<pos.size(); i++)
            { pos[i].clear();vel[i].clear(); }
        pos.clear();vel.clear();
        return boost::multi_array<float, 4>(boost::extents[0][0][0][0]);
    }
    int indz = world_rank / (nsplitx*nsplity);
    int indy = (world_rank - indz*nsplitx*nsplity) / nsplitx;
    int indx =  world_rank - indz*nsplitx*nsplity - indy*nsplitx;
    double subLx = Boxsize/nsplitx, subLy = Boxsize/nsplity, subLz = Boxsize/nsplitz;
    double subL_min[3] = { subLx*indx,     subLy*indy,     subLz*indz    }; 
    double subL_max[3] = { subLx*(indx+1), subLy*(indy+1), subLz*(indz+1) };
    if(world_size>1)
        SubDivideBox( pos, vel, subL_min, subL_max, paddingRate );

    /********* from subBox to field *********/
    DToutput output = To_Delaunay_3D( pos );
    output.field_vel = vel ;
    if( Nmesh%nsplitComm != 0 )
    {
        Nmesh += nsplitComm - Nmesh%nsplitComm;
        cerr << "warning: Grid number does match the MPI process, set to: " << Nmesh << endl;
    }
    int subNmeshx = Nmesh/nsplitx, subNmeshy = Nmesh/nsplity, subNmeshz = Nmesh/nsplitz;
    int sub_grid[3] = {subNmeshx, subNmeshy, subNmeshz};
    MeshOutput_3D<float> meshValues = InterploateGrid_3D(output, subL_min, subL_max, sub_grid );

    /********* Collect the data from all processes *********/
    int shapeSub = subNmeshx*subNmeshy*subNmeshz;
    if (world_rank != 0) 
    {
        MPI_Send(meshValues.Mesh_rho.data(), shapeSub, MPI_INT, 0, 10, MPI_COMM_WORLD);
        MPI_Send(meshValues.Mesh_vx.data(),  shapeSub, MPI_INT, 0, 11, MPI_COMM_WORLD);
        MPI_Send(meshValues.Mesh_vy.data(),  shapeSub, MPI_INT, 0, 12, MPI_COMM_WORLD);
        MPI_Send(meshValues.Mesh_vz.data(),  shapeSub, MPI_INT, 0, 13, MPI_COMM_WORLD);
        // MPI_Barrier(MPI_COMM_WORLD);
        return boost::multi_array<float, 4>(boost::extents[0][0][0][0]);
    }
    else
    {
        boost::multi_array<float, 4> arr_out(boost::extents[4][Nmesh][Nmesh][Nmesh]);
        // boost::multi_array<float, 4> subdata(boost::extents[4][subNmeshx][subNmeshy][subNmeshz]);
        // boost::multi_array<float, 3> recdata(boost::extents[subNmeshx][subNmeshy][subNmeshz]);
        vector<boost::multi_array<float, 3>> subdata(4);
        for(int itype=0; itype<4; itype++)
            subdata[itype].resize(boost::extents[subNmeshx][subNmeshy][subNmeshz]);
        subdata[0] = meshValues.Mesh_rho;
        subdata[1] = meshValues.Mesh_vx;
        subdata[2] = meshValues.Mesh_vy;
        subdata[3] = meshValues.Mesh_vz;
        // MPI_Barrier(MPI_COMM_WORLD);
        int src_rank = -1;
        for(int i2=0; i2<nsplitz; i2++)
            for(int i1=0; i1<nsplity; i1++)
                for(int i0=0; i0<nsplitx; i0++)
                {
                    src_rank += 1;
                    if(src_rank!=0)
                        for(int itype=0; itype<4; itype++)
                            MPI_Recv(subdata[itype].data(), shapeSub, MPI_INT, src_rank, 10+itype, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    int off0=i0*subNmeshx, off1=i1*subNmeshy, off2=i2*subNmeshz ;
                    for(int j0=0; j0<subNmeshx; j0++)
                        for(int j1=0; j1<subNmeshy; j1++)
                            for(int j2=0; j2<subNmeshz; j2++)
                            {
                                int k0=off0+j0, k1=off1+j1, k2=off2+j2 ;
                                for(int itype=0; itype<4; itype++)
                                    arr_out[itype][k0][k1][k2] = subdata[itype][j0][j1][j2];
                            }
                }
        return arr_out;
    }
    // MPI_Finalize();
}


