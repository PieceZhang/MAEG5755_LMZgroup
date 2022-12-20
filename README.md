# MAEG5755 Robotics - LMZgroup Project
### Robot arm manipulation and trajectory generation system
### Dec.20.2022
**Project code will be released open source at https://github.com/PieceZhang/MAEG5755_LMZgroup after final assessment** 

## Group Member
- 1155178501 LUO, Fan  
- 1155171377 MIANO, Seth Ephraim Amora  
- 1155185435 ZHANG, Yuelin  
(In alphabetical order)

## Files Description
### utilizations
- utils_equations.py: contain long equations
- utils_FK.py: methods for FK
- utils_IK.py: class for IK
- utils_trajectory_taskspace.py: class for task space trajectory planning
- utils_trajectory_jointspace.py: class for joint space trajectory planning (deprecated)
### executable
- exe_hunting_long.py: long range hunting
- exe_hunting_short.py: short range hunting
- exe_hunting_obstacle_avoid.py: short range hunting with obstacle avoid
  - (see Astar_released.m for path searching algorithm)
