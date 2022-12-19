import utils_IK as IK
import utils_trajectory_taskspace as traj

if __name__ == '__main__':
    loc_whale = [-32, 17, -5]
    loc_seagull1 = [22, 36, -4]
    loc_seagull2 = [-10, 36, -21]
    loc_end = [-20, 10, -5]

    ik = IK.IKSolverCUTER3DoFAna()
    ik.l4 = 20
    trajgen = traj.TrajQuinticContiguousTS(IKsolver=ik)
    whole_traj = trajgen.get_whole_traj([loc_end, loc_whale, loc_seagull2, loc_seagull1,
                                         loc_seagull2, loc_whale, loc_end],
                                        [0, 2, 4, 6, 7, 8, 9], [2, 4, 6])

    return_str = 'angle;'
    for trajectory in whole_traj:
        return_str += ",".join([str(i) for i in trajectory]) + ";"
    with open('long_trajectory.txt', 'w') as f:
        f.write(return_str)

    print("[INFO] Generate trajectory successfully.")
