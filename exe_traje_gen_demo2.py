from utils_trajectory_taskspace import TrajLinearTS

if __name__ == '__main__':
    # Linear trajectory (TS), which can input position sequence in any length and any format with 'fire' command
    trajgents = TrajLinearTS()
    whole_traj_ts = trajgents.get_whole_traj([[0, 0, 0], [1, 1, 1], [-2, -2, -2]], [0, 1, 2], [1, 1.5])

    exit()

    # return_str = 'angle;'
    # for trajectory in whole_traj:
    #     return_str += ",".join([str(i) for i in trajectory]) + ";"
    # with open('linear_trajectory.txt', 'w') as f:
    #     f.write(return_str)
