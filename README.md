## Dynamical Planner

``DynamicalPlanner`` is a planner determines whole-body trajectories for a humanoid robot. Check the [paper](https://ieeexplore.ieee.org/document/9196801) for a more thorough description.

### Dependencies
- [``iDynTree``](https://github.com/robotology/idyntree) (Version >= 0.11.103) It needs the  options ``IDYNTREE_COMPILES_OPTIMALCONTROL`` and ``IDYNTREE_USES_IRRLICHT`` set to ``ON``. The latter is useful for the visualization.
- [``Eigen3``](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [``matioCpp``](https://github.com/ami-iit/matio-cpp)
- [``FFmpeg``](https://ffmpeg.org/)
- [``IPOPT``](https://coin-or.github.io/Ipopt/) It is suggested to install it from source, including the HSL solvers (please check https://coin-or.github.io/Ipopt/INSTALL.html)

### Installation

The planner is compatible with Linux/macOS and Windows, but it has been mainly tested on Ubuntu 20.04.

```bash
git clone https://github.com/ami-iit/dynamical-planner
cd dynamical-planner
mkdir build && cd build
cmake ../
make
[sudo] make install
```

### Run the tests
The data for the papers has been obtained by running the [``SolverUnitTest``](https://github.com/ami-iit/dynamical-planner/blob/main/test/SolverTest.cpp) and [``SolverForComparisonUnitTest``](https://github.com/ami-iit/dynamical-planner/blob/main/test/SolverForComparisonsTest.cpp) executables. They are available after setting the ``CMake`` variable ``BUILD_TESTING`` to ``ON``.


### Cite this work

```tex
@INPROCEEDINGS{9196801,
  author={S. {Dafarra} and G. {Romualdi} and G. {Metta} and D. {Pucci}},
  booktitle={2020 IEEE International Conference on Robotics and Automation (ICRA)}, 
  title={Whole-Body Walking Generation using Contact Parametrization: A Non-Linear Trajectory Optimization Approach}, 
  year={2020},
  volume={},
  number={},
  pages={1511-1517},
  doi={10.1109/ICRA40945.2020.9196801}}
```

