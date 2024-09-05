
void MPI_Index( const int world_rank, const int world_size, int& Nprocess,
                int& indx, int& indy, int& indz )
{
    int nsplit = cbrt(world_size+1e-10);
    int indz = world_rank / pow(nsplit, 2);
    int indy = (world_rank - indz*pow(nsplit, 2)) / nsplit;
    int indx =  world_rank - indz*pow(nsplit, 2) - indy*nsplit;
    Nprocess = nsplit*nsplit*nsplit ;
}

template<typename T, typename Nd>
vector<T> FlattenToVector( boost::multi_array<T,Nd>& vecNd )
{
    vector<T> vec1d;
    if(Nd==3)
    {
        int n1 = vecNd.shape()[0], n2 = vecNd.shape()[1], n3 = vecNd.shape()[2];
        for(int i=0; i<n1; i++)
            for(int j=0; j<n2; j++)
                for(int k=0; k<n3; k++)
                    vec1d.push_back( vecNd[i][j][k] );
    }
    else if(Nd==4)
    {
        int n1 = vecNd.shape()[0], n2 = vecNd.shape()[1], n3 = vecNd.shape()[2], n4 = vecNd.shape()[3];
        for(int i=0; i<n1; i++)
            for(int j=0; j<n2; j++)
                for(int k=0; k<n3; k++)
                    for(int l=0; l<n4; l++)
                        vec1d.push_back( vecNd[i][j][k][l] );
    }
    else
        cout << "Error: Nd should be 3 or 4." << endl;
    return vec1d;
}


template<typename T, typename Nd>
void ReshapeToArray( vector<T>& vec1d, boost::multi_array<T,Nd>& vecNd, )
{
    if(Nd==3)
    {
        int n1 = vecNd.shape()[0], n2 = vecNd.shape()[1], n3 = vecNd.shape()[2];
        int start = 0;
        for(int i=0; i<n1; i++ )
            for(int j=0; j< n2;, j++ )
            {
                vec3d[i][j] = vector<T>( vec1d.begin() +start, vec1d.begin() +start+n3 );
                start += n3;
            }
    }
    else if(Nd==4)
    {
        int n1 = vecNd.shape()[0], n2 = vecNd.shape()[1], n3 = vecNd.shape()[2], n4 = vecNd.shape()[3];
        int start = 0;
        for(int i=0; i<n1; i++ )
            for(int j=0; j< n2;, j++ )
                for( int k=0; k<n3; k++ )
                {
                    vec4d[i][j][k] = vector<T>( vec1d.begin() +start, vec1d.begin() +start+n3 );
                    start += n3;
                }
    }
    else
        cout << "Error: Nd should be 3 or 4." << endl;
}