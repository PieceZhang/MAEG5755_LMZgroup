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

    def _get_traj_func(self, init_theta, end_theta, init_t, end_t):
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
            func_list = self._get_traj_func(init_theta, end_theta, timelist[num], timelist[num + 1])
            ticks_list = np.arange(timelist[num], timelist[num + 1] + 1 / self.frequency, 1 / self.frequency)
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

    def _get_traj_func(self, init_theta, end_theta, init_t, end_t):
        assert len(init_theta) == len(end_theta) == self.num_dofs
        func_list = []
        self.k = []
        self.b = []
        for init, end in zip(init_theta, end_theta):
            # self.k.append((end - init) / (end_t - init_t))
            # self.b.append(init)
            k = Symbol('k')
            b = Symbol('b')
            param = solve([k * init_t + b - init,
                           k * end_t + b - end], [k, b])
            self.k.append(float(param[k]))
            self.b.append(float(param[b]))
            func_list.append(lambda x, i: self.k[i] * x + self.b[i])
        return func_list


class TrajCubicNonContiguousTS(Trajectory):
    # Cubic trajectory
    def __init__(self, num_dofs=3, frequency=50):
        super().__init__(num_dofs=num_dofs, frequency=frequency)

    def _get_traj_func(self, init_theta, end_theta, init_t, end_t, v1=None, v2=None):
        assert len(init_theta) == len(end_theta) == self.num_dofs and type(init_theta[0]) is int
        func_list = []
        t1 = init_t
        t2 = end_t
        v1 = iter([0, 0, 0]) if v1 is None else iter(v1)
        v2 = iter([0, 0, 0]) if v2 is None else iter(v2)
        for init, end in zip(init_theta, end_theta):
            symsparam = [[Symbol('p3'), Symbol('p2'), Symbol('p1'), Symbol('p0')]]
            params = solve([symsparam[0][0] * t1 ** 3 + symsparam[0][1] * t1 ** 2 + symsparam[0][2] * t1 + symsparam[0][3] - init,
                            symsparam[0][0] * t2 ** 3 + symsparam[0][1] * t2 ** 2 + symsparam[0][2] * t2 + symsparam[0][3] - end,
                            3 * symsparam[0][0] * t1 ** 2 + 2 * symsparam[0][1] * t1 + symsparam[0][2] - v1.__next__(),
                            3 * symsparam[0][0] * t2 ** 2 + 2 * symsparam[0][1] * t2 + symsparam[0][2] - v2.__next__()],
                           symsparam[0])
            func_list.append(partial(np.polyval, p=[float(params[symsparam[0][_]]) for _ in range(4)]))
        return func_list

    def _get_values(self, t, func):
        def __dof_limit(i, x):
            if x < self.dof_min[i] or x > self.dof_max[i]:
                raise Exception("out of range!")
            else:
                return x

        return [__dof_limit(i, f(x=t)) for i, f in enumerate(func)]


class TrajCubicContiguousTS(TrajCubicNonContiguousTS):
    # Cubic trajectory
    def __init__(self, num_dofs=3, frequency=50):
        super().__init__(num_dofs=num_dofs, frequency=frequency)

    def get_whole_traj(self, thetalist: list, timelist: list, firelist=None):
        assert len(thetalist) == len(timelist)
        assert firelist is None or max(firelist) <= max(timelist)

        v_last = 0  # velocity
        value_list = [[] for i in range(self.num_dofs)]
        for num in range(len(thetalist) - 1):
            init_theta = thetalist[num]
            end_theta = thetalist[num + 1]
            if num + 1 < len(thetalist) - 1:  # intermediate point
                pass
            else:  # last point
                v = 0
            func_list = self._get_traj_func(init_theta, end_theta, timelist[num], timelist[num + 1])
            v_last = v
            ticks_list = np.arange(timelist[num], timelist[num + 1], 1 / self.frequency)
            for t in ticks_list:
                angles = self._get_values(t, func_list)
                for i in range(self.num_dofs):
                    value_list[i].append(angles[i])

        for t in firelist:
            for i in range(self.num_dofs):
                value_list[i].insert(int(t * self.frequency), 'fire')

        return value_list
