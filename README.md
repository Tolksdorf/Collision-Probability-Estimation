# Collision-Probaility-Estimation
Dependencies:  
-Eigen  
-CasaDi  

After installing CasaDi, configure the library path as: 
```
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/" >> ~/.bashrc
```
The package should be ready to be build, so go to the directory where this repo is located and run:
```
mkdir build && cd build && cmake .. && make
```
Test the probability of collision estimation with the testbench: 
```
./testbench.out
```
