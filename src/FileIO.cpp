#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <boost/multi_array.hpp>

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
    vector<vector<float>> data_float(Nsize, vector<float>(3));
    for(int i=0 ; i<Nsize ; i++)
        infile.read( reinterpret_cast<char*>(data_float[i].data()), 3* sizeof(float) );
    
    vector<vector<double>> data(Nsize, vector<double>(3));
    for(int i=0 ; i<Nsize ; i++)
        transform(data_float[i].begin(), data_float[i].end(), data[i].begin(),
                [](float f) { return static_cast<double>(f); }  );
    return data;
}


template <typename T, size_t Nd>
void SaveData( string filename, boost::multi_array<T, Nd> array )
{
    ofstream outFile(filename, ios::binary);
    outFile.write(reinterpret_cast<char*>(array.data()), array.num_elements() * sizeof(T));
}