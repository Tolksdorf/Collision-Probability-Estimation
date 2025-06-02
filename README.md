# Collision-Probability-Estimation
Tested with Ubuntu 22.04, gcc 12.1.0, and CMake 3.30.5.
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
For documentation, please see our paper on [arXiv](https://arxiv.org/abs/2505.21161). Cite it as:

```
@article{tolksdorf2025collision,
  title={Collision Probability Estimation for Optimization-based Vehicular Motion Planning},
  author={Tolksdorf, Leon and Tejada, Arturo and Birkner, Christian and van de Wouw, Nathan},
  journal={arXiv preprint arXiv:2505.21161}  [Titel anhand dieser ArXiv-ID in Citavi-Projekt übernehmen] ,
  year={2025}
}
```
