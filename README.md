# MAEG5755 Robotics - LMZgroup Project
### Robot arm manipulation and trajectory generation system
### Dec.20.2022
**written by ZHANG Yuelin**  
**Project code will be released open source at https://github.com/PieceZhang/MAEG5755_LMZgroup after final assessment** 
**See https://drive.google.com/drive/folders/1sTf8wHw-eGJKv_hgk_zIzlkV0jfu85iO?usp=share_link for test video** 

## Group Member
- LUO, Fan  
- MIANO, Seth Ephraim Amora  
- ZHANG, Yuelin  
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
