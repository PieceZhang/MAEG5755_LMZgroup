# Copyright (c) Zhang Yuelin. All Rights Reserved.
import numpy as np


class Trajectory(object):
    def __init__(self, num_joints=3, frequency=50):
        self.num_joints = num_joints
        if num_joints == 3:
            self.angle_min = [-90, -15, -180]
            self.angle_max = [90, 180, 10]
        else:
            raise ValueError("[ERROR] Undefined num_joints.")
        self.frequency = frequency  # server update frequency

    def _get_angles(self, t, func):
        """
        get angle for each joint acording to current time step
        :param t: time step
        :param func: function list of every joint,
                    each function should have 2 params, t(time) and i(num of joints)
        :return: list of lentgh 3
        """

        def __angle_limit(i, x):
            if x < self.angle_min[i] or x > self.angle_max[i]:
                raise Exception("out of range!")
            else:
                return x

        return [__angle_limit(i, f(t, i)) for i, f in enumerate(func)]

    def _get_traj_func(self, init_theta, end_theta, end_t):
        """
        get trajectory function for each joint
        Must be overriden manually
        :param init_theta: initial angle
        :param end_theta: end angle
        :param end_t: final time
        :return: list of lentgh 3
        """
        raise NotImplementedError

    def get_whole_traj(self, init_theta, end_theta, end_t):
        """
        get whole trajectory
        :param init_theta: initial angle
        :param end_theta: end angle
        :param end_t: final time
        :return: list of lentgh 3
        """
        assert len(init_theta) == len(end_theta) == self.num_joints, '[ERROR] Joints doesnt match! '

        func_list = self._get_traj_func(init_theta, end_theta, end_t)

        angle_list = [[] for i in range(self.num_joints)]
        ticks_list = np.arange(0.0, end_t, 1 / self.frequency)
        for t in ticks_list:
            angles = self._get_angles(t, func_list)
            for i in range(self.num_joints):
                angle_list[i].append(angles[i])

        return angle_list


class TrajectoryLinear(Trajectory):
    def __init__(self, num_joints=3, frequency=50):
        super().__init__(num_joints=num_joints, frequency=frequency)

    def _get_traj_func(self, init_theta, end_theta, end_t):
        func_list = []
        self.k = []
        self.b = []
        for init, end in zip(init_theta, end_theta):
            self.k.append((end - init) / end_t)
            self.b.append(init)
            func_list.append(lambda x, i: self.k[i] * x + self.b[i])
        return func_list

