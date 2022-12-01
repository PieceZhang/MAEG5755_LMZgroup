# Copyright (c) Zhang Yuelin. All Rights Reserved.
import numpy as np
from functools import partial
from sympy import Symbol, solve, nsolve, sin, cos, acos, atan, pi
import matplotlib.pyplot as plt


class _IKSolver(object):
    def __init__(self):
        """
        Basic IK solver class
        """
        pass


class _IKSolverCUTER(_IKSolver):
    def __init__(self):
        super().__init__()
        self.theta_min = [-90, 0, -140, -180, -110, -180]
        self.theta_max = [90, 180, 45, 180, 110, 180]

    def postprocessing(self, rad):
        """
        post processing, convert rad roots to angle
        :param rad: rad
        :return: angle
        """
        rad = np.array(list(map(lambda x: (x % (2 * np.pi)) / np.pi * 180, rad[0])))[None, :]  # x % (2*np.pi) ?
        # for i in range(rad.shape[1]):
        #     if rad[0, i] > self.theta_max[i]:
        #         rad[0, i] -= 360
        #         if rad[0, i] < self.theta_min[i]:
        #             print("[INFO] Warning.")
        return rad

    def __call__(self, taskspace: list):
        raise NotImplementedError


class IKSolverCUTER3DoFAna(_IKSolverCUTER):
    """
    Analytical IK
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        taskspace = np.array(taskspace).transpose()
        roots = np.ndarray((0, 3))
        rootlast = [np.pi * 0 / 180, np.pi * 8 / 180, np.pi * (-40) / 180]
        for xyz in taskspace:
            x = xyz[0]
            z = xyz[1]
            y = xyz[2]  # exchange y and z to match the coor in the simulator
            theta1 = Symbol('theta1')
            theta2 = Symbol('theta2')
            theta3 = Symbol('theta3')
            root = nsolve([(cos(theta1) * (1963 * cos(theta2 - 5361580512790693 / 36028797018963968) +
                                           2020 * cos(theta2 + theta3))) / 100 - x,
                           (sin(theta1) * (1963 * cos(theta2 - 5361580512790693 / 36028797018963968) +
                                           2020 * cos(theta2 + theta3))) / 100 - y,
                           (1963 * sin(theta2 - 5361580512790693 / 36028797018963968)) / 100 +
                           (101 * sin(theta2 + theta3)) / 5 + 509 / 50 - z],
                          [theta1, theta2, theta3], rootlast)
            rootlast = [root[0], root[1], root[2]]
            roots = np.concatenate([roots, self.postprocessing(np.array([rootlast], dtype=np.float))])
        return roots.transpose().tolist()


class IKSolverCUTER3DoFNum(_IKSolverCUTER):
    """
    Numerical IK
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        pass


class IKSolverCUTER6DoFAna(_IKSolverCUTER):
    """
    Analytical IK
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        pass


class IKSolverCUTER6DoFNum(_IKSolverCUTER):
    """
    Numerical IK
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        pass
