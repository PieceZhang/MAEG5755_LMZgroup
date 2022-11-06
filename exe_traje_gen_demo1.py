from utils_trajectory import TrajectoryLinear

if __name__ == '__main__':
    # joint initial state (0, 60, -60)
    # object1 xyz (5, 8, -35)  joint (5, 55, -100)
    # object1 xyz (0, 50 -20)  joint (0, 105, -60)

    trajgen = TrajectoryLinear()
    whole_traj = trajgen.get_whole_traj([-90, 100, -100], [5, 55, -100], 3)
    # whole_traj = trajgen.get_whole_traj([-90, 100, -100], [0, 105, -60], 3)

    whole_traj[0] = whole_traj[0] + ['fire'] + list(reversed(whole_traj[0])) + ['fire']
    whole_traj[1] = whole_traj[1] + ['fire'] + list(reversed(whole_traj[1])) + ['fire']
    whole_traj[2] = whole_traj[2] + ['fire'] + list(reversed(whole_traj[2])) + ['fire']

    return_str = 'angle;'
    for trajectory in whole_traj:
        return_str += ",".join([str(i) for i in trajectory]) + ";"
    with open('linear_trajectory.txt', 'w') as f:
        f.write(return_str)
