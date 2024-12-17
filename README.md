# CWF: Consolidating Weak Features in High-quality Mesh Simplification

Code of CWF: Consolidating Weak Features in High-quality Mesh Simplification.

Project Page: https://ruixu.me/html/CWF/index.html 

Paper link: https://arxiv.org/abs/2404.15661


Note: If you have trouble configuring, please try our [original version](https://github.com/Xrvitd/CWF/tree/c2d6a7fe35bac171063eb77e66df9dc9eca8a2d8), as well as this version of [vcpkg](https://www.dropbox.com/scl/fi/x79za5w28sg4yzq53uh7e/vcpkg.zip?rlkey=re26jx3rfb6nw2gv776w6b0k9&dl=0). 

Please cite our work:
```
@article{xu2024cwf,
      title={CWF: Consolidating Weak Features in High-quality Mesh Simplification}, 
      author={Xu, Rui and Liu, Longdu and Wang, Ningna and Chen, Shuangmin and Xin, Shiqing and Guo, Xiaohu and Zhong, Zichun and Komura, Taku and Wang, Wenping and Tu, Changhe},
      journal={ACM Transactions on Graphics (TOG)},
      publisher={ACM New York, NY, USA},
      year={2024},
      address = {New York, NY, USA},
      volume = {43},
      number = {4},
      issn = {0730-0301},
      url = {https://doi.org/10.1145/3658159},
      doi = {10.1145/3658159},
}
```

Currently we have only tested our code on 64-bit windows systems and Visual Studio 2022 Professional.

### Dependence

- CGAL
- Eigen3
- Boost
- Libigl

### Please using vcpkg to install dependent libraries!!!

#### Important: Please use  "git clone" to install vcpkg, otherwise you may get errors in cgal installation.

- .\vcpkg install boost:x64-windows
- .\vcpkg install cgal:x64-windows
  
  ​	use "git pull" if you get errors with the "gmp" library.
- .\vcpkg install Eigen3:x64-windows
- .\vcpkg install libigl:x64-windows
- .\vcpkg integrate install

### MSVC on Windows

```
Download this project: CWF
```

Open Cmake-GUI

```
Where is the source code: CWF

Where to build the binaries: CFW/build
```

note: check the location of dependencies and install. It is recommended to use vcpkg to add dependencies.

Configure->Generate->Open Project

Then, build the project:

ALL_BUILD -> Build
Turn Debug to Release -> ALL_BUILD -> Build


## Test

The example is in 'MAIN'.

All the files is in 'CWF\data'.

The output files is in 'CWF\data\LBFGSOUT'.

We put the result of our operation in 'CWF\data\LBFGSOUT\DemoOutput', you can use it for comparison to know whether the program is running correctly.

## IMPORTANT NOTE:

Make sure you are using C++ 14, not C++ 17 or newer versions. 

This code is not optimized for speed, but for clarity.

The default number of Openmp parallel threads is 30, set according to an AMD Ryzen 5950x CPU, please set different number of threads according to the CPU you use to get the best running effect.



## Testing Platform

- Windows 10
- Visual Studio 2022 Professional
- AMD Ryzen 5950X
- 64GB Memory

Considering that most computers may not have this configuration, this commit does not support acceleration.

In order to allow you to view the optimization process in more detail,
we have only set the optimization stop at 50 iterations, you can manually stop the optimization, and view all iteration results in the data\LBFGSOUT folder.


## License
CWF is under [AGPL-3.0](https://www.gnu.org/licenses/agpl-3.0.en.html), so any downstream solution and products (including cloud services) that include CWF code inside it should be open-sourced to comply with the AGPL conditions. For learning purposes only and not for commercial use. If you want to use it for commercial purposes, please contact us first.

