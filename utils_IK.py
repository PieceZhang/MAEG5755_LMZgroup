# Copyright (c) Zhang Yuelin. All Rights Reserved.
import math
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
        for i in range(rad.shape[1]):
            if rad[0, i] > self.theta_max[i]:
                rad[0, i] -= 360
                if rad[0, i] < self.theta_min[i]:
                    print("[INFO] Warning.")
            elif rad[0, i] < self.theta_min[i]:
                rad[0, i] += 360
                if rad[0, i] > self.theta_max[i]:
                    print("[INFO] Warning.")
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
        rootlast = [np.pi * 0 / 180, np.pi * (-0) / 180]
        for xyz in taskspace:
            x = xyz[0]
            z = xyz[1]
            y = xyz[2]  # exchange y and z to match the coor in the simulator
            # theta1 = Symbol('theta1')
            theta1 = atan(y / x)
            theta2 = Symbol('theta2')
            theta3 = Symbol('theta3')
            # A = (19.63 * cos(theta2 - 93 / 625) + 20.2 * cos(theta2 + theta3))
            # B = (19.63 * sin(theta2 - 93 / 625) + 20.2 * sin(theta2 + theta3))
            A = -x * sin(theta1) + y * cos(theta1)
            B = z - 10.18
            root = nsolve([acos((A ** 2 + B ** 2 - 19.63 ** 2 - 20.2 ** 2) / (2 * 19.63 * 20.2)) - theta3,
                           atan((B * (19.63 + 20.2 * cos(theta3 + 93 / 625)) - A * 20.2 * sin(theta3 + 93 / 62)) /
                                (A * (19.63 + 20.2 * cos(theta3 + 93 / 625)) + B * 20.2 * sin(theta3 + 93 / 62))) - theta2],
                          [theta2, theta3], rootlast)
            rootlast = [root[0], root[1]]
            if len(rootlast[0].args) == 2:
                print(len(rootlast[0].args))
                roots = np.concatenate([roots, self.postprocessing(np.array([[theta1, rootlast[0].args[0], rootlast[1].args[0]]], dtype=np.float))])
            else:
                roots = np.concatenate([roots, self.postprocessing(np.array([[theta1, root[0], root[1]]], dtype=np.float))])
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
