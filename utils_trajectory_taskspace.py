# Copyright (c) Zhang Yuelin. All Rights Reserved.
import numpy as np
from functools import partial
from sympy import Symbol, solve
import matplotlib.pyplot as plt


class _Trajectory(object):
    def __init__(self, IKsolver=None, num_dofs=3, frequency=50):
        """
        Basic trajectory class
        :param num_dofs: number of DoFs (q, task space)
        :param frequency: simulation freq
        """
        self.IKsolver = IKsolver
        self.num_dofs = num_dofs
        if num_dofs == 3:
            self.dof_min = [-60, -60, -60]  # TBD
            self.dof_max = [60, 60, 60]
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

    def get_whole_traj(self, thetalist: list, timelist: list, firelist=None, visual=False):
        """
        get whole trajectory
        :param visual:
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

        if visual:
            for i in range(self.num_dofs):
                plt.plot(np.arange(0, len(value_list[i]), 1), value_list[i])
            plt.show()

        if self.IKsolver is not None:
            value_list = self.IKsolver(value_list)

        for t in firelist:
            for i in range(self.num_dofs):
                value_list[i].insert(int(t * self.frequency), 'fire')

        return value_list


class TrajLinearTS(_Trajectory):
    """
    Linear trajectory
    """
    def _get_traj_func(self, init_theta, end_theta, init_t, end_t):
        assert len(init_theta) == len(end_theta) == self.num_dofs
        func_list = []
        self.k = []
        self.b = []
        for init, end in zip(init_theta, end_theta):
            k = Symbol('k')
            b = Symbol('b')
            param = solve([k * init_t + b - init,
                           k * end_t + b - end], [k, b])
            self.k.append(float(param[k]))
            self.b.append(float(param[b]))
            func_list.append(lambda x, i: self.k[i] * x + self.b[i])
        return func_list


class TrajCubicNonContiguousTS(_Trajectory):
    """
    Cubic trajectory
    """
    def _get_traj_func(self, init_theta, end_theta, init_t, end_t, v1=None, v2=None):
        assert len(init_theta) == len(end_theta) == self.num_dofs and type(init_theta[0]) is int
        func_list = []
        t1 = init_t
        t2 = end_t
        v1 = iter([0 for _ in range(self.num_dofs)]) if v1 is None else iter(v1)
        v2 = iter([0 for _ in range(self.num_dofs)]) if v2 is None else iter(v2)
        for init, end in zip(init_theta, end_theta):
            symsparam = [Symbol('p3'), Symbol('p2'), Symbol('p1'), Symbol('p0')]
            params = solve([symsparam[0] * t1 ** 3 + symsparam[1] * t1 ** 2 + symsparam[2] * t1 + symsparam[3] - init,
                            symsparam[0] * t2 ** 3 + symsparam[1] * t2 ** 2 + symsparam[2] * t2 + symsparam[3] - end,
                            3 * symsparam[0] * t1 ** 2 + 2 * symsparam[1] * t1 + symsparam[2] - v1.__next__(),
                            3 * symsparam[0] * t2 ** 2 + 2 * symsparam[1] * t2 + symsparam[2] - v2.__next__()],
                           symsparam)
            func_list.append(partial(np.polyval, p=[float(params[symsparam[_]]) for _ in range(len(symsparam))]))
        return func_list

    def _get_values(self, t, func):
        def __dof_limit(i, x):
            if x < self.dof_min[i] or x > self.dof_max[i]:
                raise Exception("out of range!")
            else:
                return x

        return [__dof_limit(i, f(x=t)) for i, f in enumerate(func)]


class TrajCubicContiguousTS(TrajCubicNonContiguousTS):
    """
    Cubic trajectory
    """
    def get_whole_traj(self, thetalist: list, timelist: list, firelist=None, visual=False):
        assert len(thetalist) == len(timelist)
        assert firelist is None or max(firelist) <= max(timelist)

        v_last = [0 for _ in range(self.num_dofs)]  # velocity
        v = [0 for _ in range(self.num_dofs)]
        value_list = [[] for i in range(self.num_dofs)]
        for num in range(len(thetalist) - 1):
            init_theta = thetalist[num]
            end_theta = thetalist[num + 1]
            if num + 1 < len(thetalist) - 1:  # intermediate point
                for i in range(self.num_dofs):
                    if np.sign(end_theta[i] - init_theta[i]) == np.sign(thetalist[num + 2][i] - end_theta[i]):  # if in the same direction
                        v[i] = (thetalist[num + 2][i] - init_theta[i]) / (timelist[num + 2] - timelist[num])  # midpoint velocity
                    else:
                        v[i] = 0
            else:  # last point
                v = [0 for _ in range(self.num_dofs)]
            func_list = self._get_traj_func(init_theta, end_theta, timelist[num], timelist[num + 1], v_last, v)
            v_last = v.copy()
            ticks_list = np.arange(timelist[num], timelist[num + 1] + 1 / self.frequency, 1 / self.frequency)
            for t in ticks_list:
                angles = self._get_values(t, func_list)
                for i in range(self.num_dofs):
                    value_list[i].append(angles[i])
            for i in range(self.num_dofs):
                value_list[i][-1] = np.round(value_list[i][-1])  # eliminate residue

        if visual:
            for i in range(self.num_dofs):
                plt.plot(np.arange(0, len(value_list[i]), 1), value_list[i])
            plt.show()

        if self.IKsolver is not None:
            value_list = self.IKsolver(value_list)

        for t in firelist:
            for i in range(self.num_dofs):
                value_list[i].insert(int(t * self.frequency), 'fire')

        return value_list


class TrajQuinticContiguousTS(TrajCubicContiguousTS):
    """
    Quintic trajectory
    """
    def _get_traj_func(self, init_theta, end_theta, init_t, end_t, v1=None, v2=None, a1=None, a2=None):
        assert len(init_theta) == len(end_theta) == self.num_dofs and type(init_theta[0]) is int
        func_list = []
        t1 = init_t
        t2 = end_t
        v1 = iter([0 for _ in range(self.num_dofs)]) if v1 is None else iter(v1)
        v2 = iter([0 for _ in range(self.num_dofs)]) if v2 is None else iter(v2)
        a1 = iter([0 for _ in range(self.num_dofs)]) if a1 is None else iter(a1)
        a2 = iter([0 for _ in range(self.num_dofs)]) if a2 is None else iter(a2)
        for init, end in zip(init_theta, end_theta):
            symsparam = [Symbol('p5'), Symbol('p4'), Symbol('p3'), Symbol('p2'), Symbol('p1'), Symbol('p0')]
            params = solve([symsparam[0] * t1 ** 5 + symsparam[1] * t1 ** 4 + symsparam[2] * t1 ** 3 +
                            symsparam[3] * t1 ** 2 + symsparam[4] * t1 + symsparam[5] - init,
                            symsparam[0] * t2 ** 5 + symsparam[1] * t2 ** 4 + symsparam[2] * t2 ** 3 +
                            symsparam[3] * t2 ** 2 + symsparam[4] * t2 + symsparam[5] - end,
                            5 * symsparam[0] * t1 ** 4 + 4 * symsparam[1] * t1 ** 3 + 3 * symsparam[2] * t1 ** 2 +
                            2 * symsparam[3] * t1 + symsparam[4] - v1.__next__(),
                            5 * symsparam[0] * t2 ** 4 + 4 * symsparam[1] * t2 ** 3 + 3 * symsparam[2] * t2 ** 2 +
                            2 * symsparam[3] * t2 + symsparam[4] - v2.__next__(),
                            20 * symsparam[0] * t1 ** 3 + 12 * symsparam[1] * t1 ** 2 + 6 * symsparam[2] * t1 + 2 * symsparam[3] - a1.__next__(),
                            20 * symsparam[0] * t2 ** 3 + 12 * symsparam[1] * t2 ** 2 + 6 * symsparam[2] * t2 + 2 * symsparam[3] - a2.__next__()],
                           symsparam)
            func_list.append(partial(np.polyval, p=[float(params[symsparam[_]]) for _ in range(len(symsparam))]))
        return func_list


"""
# Define a trajectory generator on your own

class TrajDIY(_Trajectory):

    def __init__(self, IKsolver=None, num_dofs=3, frequency=50):
        super().__init__(IKsolver=None, num_dofs=num_dofs, frequency=frequency)

    def _get_traj_func(self, init_theta, end_theta, init_t, end_t):
        assert len(init_theta) == len(end_theta) == self.num_dofs
        ...

"""