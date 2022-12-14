from utils_trajectory_taskspace import TrajLinearTS, TrajCubicNonContiguousTS, TrajCubicContiguousTS, TrajQuinticContiguousTS

if __name__ == '__main__':
    # Linear trajectory (TS), which can input position sequence in any length and any format with 'fire' command
    trajgents = TrajLinearTS()
    whole_traj_ts_l = trajgents.get_whole_traj([[2, 3, 4], [1, 1, 1], [-2, -2, -2], [5, 5, 5]], [0, 1, 2, 2.5], [1, 1.5])

    # Cubic trajectory (TS), which can input position sequence in any length and any format with 'fire' command
    trajgents = TrajCubicNonContiguousTS()
    whole_traj_ts_nc = trajgents.get_whole_traj([[2, 3, 4], [-1, -1, -1], [1, 1, 1], [5, 5, 5], [6, 6, 6], [9, 8, 7]],
                                                [0, 1, 2, 2.5, 3, 5], [1, 1.5], True)

    # Cubic trajectory (TS), which can input position sequence in any length and any format with 'fire' command
    trajgents = TrajCubicContiguousTS()
    whole_traj_ts_c = trajgents.get_whole_traj([[2, 3, 4], [-1, -1, -1], [1, 1, 1], [5, 5, 5], [6, 6, 6], [9, 8, 7]],
                                               [0, 1, 2, 2.5, 3, 5], [1, 1.5], True)

    # Quintic trajectory (TS), which can input position sequence in any length and any format with 'fire' command
    trajgents = TrajQuinticContiguousTS()
    whole_traj_ts_q = trajgents.get_whole_traj([[2, 3, 4], [-1, -1, -1], [1, 1, 1], [5, 5, 5], [6, 6, 6], [9, 8, 7]],
                                               [0, 1, 2, 2.5, 3, 5], [1, 1.5], True)

    exit()

    # return_str = 'angle;'
    # for trajectory in whole_traj:
    #     return_str += ",".join([str(i) for i in trajectory]) + ";"
    # with open('linear_trajectory.txt', 'w') as f:
    #     f.write(return_str)
