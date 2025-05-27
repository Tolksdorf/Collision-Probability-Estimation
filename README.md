# Collision-Probability-Estimation
Dependencies:  
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)

  ```
  sudo apt install -y libeigen3-dev
  ```
- CasADi, built from source, see [here](https://github.com/casadi/casadi/wiki/InstallationLinux).  

After installing CasaDi, configure the library path as: 
```
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/" >> $HOME/.bashrc && \
source $HOME/.bashrc
```
The package is ready to be built, so go to the directory where this repo is located and run:
```
mkdir build && cd build && cmake .. && make
```
Test the probability of collision estimation with the testbench: 
```
./testbench.out
```
For documentation, please see our paper on arXiv.
