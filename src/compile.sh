# ------------------------------------------------------------------
# load the necessary modules
module load compiler/gcc-11.3.0
module load boost/1_73_0
module load compiler/cmake-3.25.1

# ------------------------------------------------------------------
## The cmake should be in the same environment when installing the CGAL-package, since it requires the same gcc/g++ compiler.
## $module load anaconda/anaconda-mamba
## It does not properly load the conda environment in the subshells.
source /opt/conda/conda-4.12.0/etc/profile.d/conda.sh
conda activate shuren_env3.9


#g++ -o out main.cpp -I/home/shuren/Code/CGAL-5.6.1/MyInstall/include -std=c++17
rm -rf ./build/*
mkdir build
cd build
cmake ..
cmake --build . --config Release # --target check
if [ $? -eq 0 ]; then
    cd ..
else
    echo "Compilation failed"
    exit 1
fi
