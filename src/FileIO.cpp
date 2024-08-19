#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;



inline bool FileExists(string& path)
{
    ifstream f( path, ios::binary );
    if( !f.good() )
    {
        cout << "File: " << path << "  Does NOT Exist." << endl;
        exit(0);
    }
    return true;
}


vector<vector<double>> ReadData( string& path, int Nsize )
{
    ifstream infile( path , ios::binary );
    vector<vector<float>> data_float(3);
    for( int j=0; j<3; data_float[j++].resize( Nsize ) );
    infile.read( reinterpret_cast<char*>(data_float[0].data()), Nsize* sizeof(float) );
    infile.read( reinterpret_cast<char*>(data_float[1].data()), Nsize* sizeof(float) );
    infile.read( reinterpret_cast<char*>(data_float[2].data()), Nsize* sizeof(float) );
    
    vector<vector<double>> data(3);
    for( int j=0; j<3; data[j++].resize( Nsize ) );
    for(int i=0 ; i<3 ; i++)
        transform(data_float[i].begin(), data_float[i].end(), data[i].begin(),
                [](float f) { return static_cast<double>(f); }  );
    return data;
}

vector<vector<double>> ReadData2D( string& path, int Nsize )
{
    ifstream infile( path , ios::binary );
    vector<vector<float>> data_float(2);
    for( int j=0; j<2; data_float[j++].resize( Nsize ) );
    infile.read( reinterpret_cast<char*>(data_float[0].data()), Nsize* sizeof(float) );
    infile.read( reinterpret_cast<char*>(data_float[1].data()), Nsize* sizeof(float) );
    
    vector<vector<double>> data(2);
    for( int j=0; j<2; data[j++].resize( Nsize ) );
    for(int i=0 ; i<2 ; i++)
        transform(data_float[i].begin(), data_float[i].end(), data[i].begin(),
                [](float f) { return static_cast<double>(f); }  );
    return data;
}



void SaveData( string filename, vector<vector<vector<double>>> arr3 )
{
    ofstream outFile(filename, ios::binary);
    for( auto &arr2 : arr3)
        for( auto &arr1 : arr2)
            outFile.write(reinterpret_cast<const char*>(arr1.data()), arr1.size()*sizeof(double));
    outFile.close();
}

void SaveData( string filename, vector<vector<double>> arr2 )
{
    ofstream outFile(filename, ios::binary);
    for( auto &arr1 : arr2)
        outFile.write(reinterpret_cast<const char*>(arr1.data()), arr1.size()*sizeof(double));
    outFile.close();
}