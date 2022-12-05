import utils_IK as IK
import utils_trajectory_taskspace as traj

if __name__ == '__main__':
    loc_pigeon = [5, 8, -35]
    loc_Penguin = [-12.585, 10.906, -45.686]
    loc_box3 = [-21.496, -1.429, -24.649]
    loc_box2 = [-10.055, 4.032, -30.630]
    loc_box1 = [-11.709, -0.302, -44.058]
    loc_seagull1 = [0.000, 40.000, -20.000]
    loc_seagull2 = [-20.116, 47.255, -6.148]
    loc_end = [-25, 20, -20]

    trajgen = traj.TrajLinearTS(IK.IKSolverCUTER3DoFAna())
    # whole_traj = trajgen.get_whole_traj([[-4, 3, -32], [-25, 20, -10], [-4, 3, -32]], [0, 1, 2], [1, 2])  # -10 30 -80
    # whole_traj = trajgen.get_whole_traj([loc_end, loc_pigeon, loc_end, loc_pigeon], [0, 1, 2, 3], [1, 2])
    # whole_traj = trajgen.get_whole_traj([loc_end, loc_Penguin, loc_end, loc_Penguin], [0, 1, 2, 3], [1, 2])
    # whole_traj = trajgen.get_whole_traj([loc_end, loc_box1, loc_end, loc_box1], [0, 1, 2, 3], [1, 2])
    whole_traj = trajgen.get_whole_traj([loc_end, loc_box2, loc_end, loc_box2], [0, 1, 2, 3], [1, 2])
    # whole_traj = trajgen.get_whole_traj([loc_end, loc_box3, loc_end, loc_box3], [0, 1, 2, 3], [1, 2])
    # whole_traj = trajgen.get_whole_traj([loc_end, loc_seagull1, loc_end, loc_seagull1], [0, 1, 2, 3], [1, 2])
    # whole_traj = trajgen.get_whole_traj([loc_end, loc_seagull2, loc_end, loc_seagull2], [0, 1, 2, 3], [1, 2])

    return_str = 'angle;'
    for trajectory in whole_traj:
        return_str += ",".join([str(i) for i in trajectory]) + ";"
    with open('linear_trajectory.txt', 'w') as f:
        f.write(return_str)
