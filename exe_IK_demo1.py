import utils_IK as IK
import utils_trajectory_taskspace as traj

if __name__ == '__main__':
    # q   (-20, 90, -90) ->
    # xyz (-7, 30, -22) -> (-2.726, 3.528, -15.606)

    trajgen = traj.TrajLinearTS(IK.IKSolverCUTER3DoFAna())
    whole_traj = trajgen.get_whole_traj([[-7, 30, -22], [0, 5.000, -26.523], [-7, 30, -22]], [0, 2, 4], [2, 4])

    return_str = 'angle;'
    for trajectory in whole_traj:
        return_str += ",".join([str(i) for i in trajectory]) + ";"
    with open('linear_trajectory.txt', 'w') as f:
        f.write(return_str)
