import utils_IK as IK
import utils_trajectory_taskspace as traj

if __name__ == '__main__':

    trajgen = traj.TrajLinearTS(IK.IKSolverCUTER3DoFAna())
    whole_traj = trajgen.get_whole_traj([[5, 8, -35], [-25, 20, -10], [5, 8, -35]], [0, 1, 2], [0, 1])
    # whole_traj = trajgen.get_whole_traj([[-4, 3, -32], [-25, 20, -10], [-4, 3, -32]], [0, 1, 2], [1, 2])  # -10 30 80

    return_str = 'angle;'
    for trajectory in whole_traj:
        return_str += ",".join([str(i) for i in trajectory]) + ";"
    with open('linear_trajectory.txt', 'w') as f:
        f.write(return_str)
