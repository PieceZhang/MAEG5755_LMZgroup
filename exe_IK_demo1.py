import utils_IK as IK
import utils_trajectory_taskspace as traj

if __name__ == '__main__':
    # (-23, 30, 0) -> (-2.726, 3.528, -15.606)

    trajgen = traj.TrajLinearTS(IK.IKSolverCUTER3DoFAna())
    whole_traj = trajgen.get_whole_traj([[-23, 30, 0], [-2.726, 3.528, -15.606], [-23, 30, 0]], [0, 2, 4], [2, 4])

    return_str = 'angle;'
    for trajectory in whole_traj:
        return_str += ",".join([str(i) for i in trajectory]) + ";"
    with open('linear_trajectory.txt', 'w') as f:
        f.write(return_str)
