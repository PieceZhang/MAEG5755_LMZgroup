# Copyright (c) Zhang Yuelin. All Rights Reserved.
import numpy as np
from functools import partial
from sympy import Symbol, solve


class Trajectory(object):
    def __init__(self, num_dofs=3, frequency=50):
        self.num_dofs = num_dofs
        if num_dofs == 3:
            self.dof_min = [-10, -10, -10]  # TODO
            self.dof_max = [10, 10, 10]
        else:
            raise ValueError("[ERROR] Undefined num_dofs.")
        self.frequency = frequency  # server update frequency

    def _get_values(self, t, func):
        """
        get value for each dof acording to current time step
        :param t: time step
        :param func: function list of every joint,
                    each function should have 2 params, t(time) and i(num of dofs)
        :return: list
        """

        def __dof_limit(i, x):
            if x < self.dof_min[i] or x > self.dof_max[i]:
                raise Exception("out of range!")
            else:
                return x

        return [__dof_limit(i, f(t, i)) for i, f in enumerate(func)]

    def _get_traj_func(self, init_theta, end_theta, end_t):
        """
        get trajectory function for each dof
        Must be overriden manually
        :param init_theta: initial angle
        :param end_theta: end angle
        :param end_t: final time
        :return: list of lentgh 3
        """
        raise NotImplementedError

    def get_whole_traj(self, thetalist: list, timelist: list, firelist=None):
        """
        get whole trajectory
        :param firelist: fire command list
        :param thetalist: list of dof
        :param timelist: time list
        :return: list
        """
        assert len(thetalist) == len(timelist)
        assert firelist is None or max(firelist) <= max(timelist)

        value_list = [[] for i in range(self.num_dofs)]
        for num in range(len(thetalist) - 1):
            init_theta = thetalist[num]
            end_theta = thetalist[num + 1]
            end_t = timelist[num + 1] - timelist[num]  # time duation of current step
            func_list = self._get_traj_func(init_theta, end_theta, end_t)
            ticks_list = np.arange(0.0, end_t, 1 / self.frequency)
            for t in ticks_list:
                angles = self._get_values(t, func_list)
                for i in range(self.num_dofs):
                    value_list[i].append(angles[i])

        for t in firelist:
            for i in range(self.num_dofs):
                value_list[i].insert(int(t * self.frequency), 'fire')

        return value_list


class TrajLinearTS(Trajectory):
    # Linear trajectory
    def __init__(self, num_dofs=3, frequency=50):
        super().__init__(num_dofs=num_dofs, frequency=frequency)

    def _get_traj_func(self, init_theta, end_theta, end_t):
        assert len(init_theta) == len(end_theta) == self.num_dofs
        func_list = []
        self.k = []
        self.b = []
        for init, end in zip(init_theta, end_theta):
            self.k.append((end - init) / end_t)
            self.b.append(init)
            func_list.append(lambda x, i: self.k[i] * x + self.b[i])
        return func_list


# class TrajCubicNonContiguousTS(Trajectory):
#     # Cubic trajectory
#     def __init__(self, num_dofs=3, frequency=50):
#         super().__init__(num_dofs=num_dofs, frequency=frequency)
#
#     def _get_traj_func(self, init_theta, end_theta, end_t):
#         assert len(init_theta) == len(end_theta) == self.num_dofs
#         func_list = []
#         if type(init_theta[0]) is int:
#             symsparam = [[Symbol('p%d3' % i), Symbol('p%d2' % i), Symbol('p%d1' % i), Symbol('p%d0' % i)] for i in [0]]
#         else:
#             raise TypeError()
#
#         for init, end in zip(init_theta, end_theta):
#             solve([])
#             func_list.append()
#         return func_list
#
#
# class TrajCubicContiguousTS(Trajectory):
#     # Cubic trajectory
#     def __init__(self, num_dofs=3, frequency=50):
#         super().__init__(num_dofs=num_dofs, frequency=frequency)
#
#     def _get_traj_func(self, init_theta, end_theta, end_t):
#         assert len(init_theta) == len(end_theta) == self.num_dofs
#         func_list = []
#         if type(init_theta[0]) is list:
#             symsparam = [[Symbol('p%d3' % i), Symbol('p%d2' % i), Symbol('p%d1' % i), Symbol('p%d0' % i)] for i in range(len(init_theta[0]))]
#         else:
#             raise TypeError()
#
#         for init, end in zip(init_theta, end_theta):
#             func_list.append()
#         return func_list

