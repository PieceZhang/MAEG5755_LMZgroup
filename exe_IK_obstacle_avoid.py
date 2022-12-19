import utils_IK as IK
import utils_trajectory_taskspace as traj

if __name__ == '__main__':
    loc_Penguin_right = [-6, 6, -35]
    path_Penguin_right = [[-8, 30, -8],
                          [-2, 26.5, -11.5],
                          [-0.50, 23, -15],
                          [0.5, 19.5, -18.5],
                          [2.50, 16.0, -22],
                          [2.50, 16.0, -25.5],
                          [-1.00, 12.5, -29.0],
                          [-4.50, 9.00, -32.5]]

    loc_Penguin_left = [24, 5, -26]
    path_Penguin_left = [[-8, 30, -8],
                         [-1.00000000000000, 23, -11.5000000000000],
                         [6.00000000000000, 16.0000000000000, -18.5000000000000],
                         [13.0000000000000, 12.5000000000000, -22],
                         [20, 5.50000000000000, -25.5000000000000]]

    # loc_end = [6, 25, - 32]

    trajgen = traj.TrajQuinticContiguousTS(IKsolver=IK.IKSolverCUTER6DoFNumxyz())
    # whole_traj = trajgen.get_whole_traj(path_Penguin_right + [loc_Penguin_right] + list(reversed(path_Penguin_right)),
    #                                     [0.2 * _ for _ in range(len(path_Penguin_right) * 2 + 1)],
    #                                     [0.2 * (len(path_Penguin_right)), 0.2 * (len(path_Penguin_right) * 2)])
    # whole_traj = trajgen.get_whole_traj(path_Penguin_right + [loc_Penguin_right, loc_end],
    #                                     [0.2 * _ for _ in range(len(path_Penguin_right) + 1)]
    #                                     + [0.2 * (len(path_Penguin_right)) + 1],
    #                                     [0.2 * (len(path_Penguin_right)), 0.2 * (len(path_Penguin_right)) + 1])
    whole_traj = trajgen.get_whole_traj(path_Penguin_left + [loc_Penguin_left] + list(reversed(path_Penguin_left)),
                                        [0.2 * _ for _ in range(len(path_Penguin_left) * 2 + 1)],
                                        [0.2 * (len(path_Penguin_left))])

    return_str = 'angle;'
    for trajectory in whole_traj:
        return_str += ",".join([str(i) for i in trajectory]) + ";"
    with open('linear_trajectory.txt', 'w') as f:
        f.write(return_str)

    print("[INFO] Generate trajectory successfully.")
