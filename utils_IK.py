# Copyright (c) Zhang Yuelin. All Rights Reserved.
import math

import numpy as np
from functools import partial
from numpy import sin, cos, arccos, arcsin, arctan2, pi, sqrt
# from sympy import Symbol, solve, nsolve, sin, cos, acos, atan, pi
import matplotlib.pyplot as plt


def rad2deg(x):
    return x / pi * 180


def deg2rad(x):
    return x / 180 * pi


class _IKSolverCUTER(object):
    def __init__(self):
        """
        Basic IK solver class
        """
        pass


class _IKSolverCUTER3DoF(_IKSolverCUTER):
    def __init__(self):
        super().__init__()
        self.theta_min = [-90, -15, -140]
        self.theta_max = [90, 180, 45]
        self.l1 = 10.18
        self.l2 = 19.41
        self.l3 = 2.91
        self.l4 = 20.2

    def postprocessing(self, rad):
        """
        post processing, convert rad roots to angle
        :param rad: rad
        :return: angle
        """
        deg = list(map(lambda x: rad2deg(x), rad))
        deg[1] = -deg[1] + 240  # theta2 offset
        constrain = [[[self.theta_min[0] < deg[0][0, 0] < self.theta_max[0],
                       self.theta_min[0] < deg[0][0, 1] < self.theta_max[0]],
                      [self.theta_min[0] < deg[0][1, 0] < self.theta_max[0],
                       self.theta_min[0] < deg[0][1, 1] < self.theta_max[0]]],
                     [[self.theta_min[1] < deg[1][0, 0] < self.theta_max[1],
                       self.theta_min[1] < deg[1][0, 1] < self.theta_max[1]],
                      [self.theta_min[1] < deg[1][1, 0] < self.theta_max[1],
                       self.theta_min[1] < deg[1][1, 1] < self.theta_max[1]]],
                     [[self.theta_min[2] < deg[2][0] < self.theta_max[2],
                       self.theta_min[2] < deg[2][1] < self.theta_max[2]]]]
        for theta3 in range(2):
            # theta3
            if not constrain[2][0][theta3]:
                continue
            # theta2
            try:
                theta2 = constrain[1][theta3].index(True)
            except ValueError:
                continue
            # theta1
            if not constrain[0][theta3][theta2]:
                continue
                # if not constrain[0][theta3][int(not theta2)]:
                #     continue
                # else:
                #     theta1 = int(not theta2)
                #     break
            else:
                theta1 = theta2
                break
        if not constrain[2][0][theta3] or not ('theta2' in locals().keys()) or not ('theta1' in locals().keys()):
            raise ValueError('[IKSolver] Location cannot reach!')
        return np.array([[deg[0][theta3, theta1], deg[1][theta3, theta2], deg[2][theta3]]], dtype=float)

    def __call__(self, taskspace: list):
        raise NotImplementedError


class IKSolverCUTER3DoFAna(_IKSolverCUTER3DoF):
    """
    Analytical IK
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        # redefine variables name
        l1 = self.l1
        l2 = sqrt(self.l2 ** 2 + self.l3 ** 2)
        l3 = self.l4
        taskspace = np.array(taskspace).transpose()
        qt = np.ndarray((0, 3))
        for xyz in taskspace:
            x = xyz[0]
            z = xyz[1]
            y = xyz[2]  # exchange y and z to match the coor in the simulator
            # solving for theta3
            theta3 = arccos((x ** 2 + y ** 2 + (z - l1) ** 2 - l2 ** 2 - l3 ** 2) / (2 * l2 * l3))
            theta3 = np.array([theta3, -theta3])
            # solving for theta2
            a = l2 + l3 * cos(theta3)
            b = l3 * cos(theta3)
            alpha = arctan2(b / sqrt(a ** 2 + b ** 2), a / sqrt(a ** 2 + b ** 2))
            asi = arcsin((z - l1) / sqrt(a ** 2 + b ** 2))
            theta2 = np.array([[asi[0] - alpha[0], pi - (asi[0] - alpha[0])],
                               [asi[1] - alpha[1], pi - (asi[1] - alpha[1])]])
            # solving for theta1
            sin1_0 = -x / (l2 * cos(theta2[0]) + l3 * cos(theta2[0] + theta3[0]))
            cos1_0 = y / (l2 * cos(theta2[0]) + l3 * cos(theta2[0] + theta3[0]))
            sin1_1 = -x / (l2 * cos(theta2[1]) + l3 * cos(theta2[1] + theta3[1]))
            cos1_1 = y / (l2 * cos(theta2[1]) + l3 * cos(theta2[1] + theta3[1]))
            theta1 = np.array([arctan2(sin1_1, cos1_1), arctan2(sin1_0, cos1_0)])  # 1, 0 ?
            # cat
            theta3 -= 0.1488
            theta2 += 0.1488
            # theta2 += np.array([[0.1488, -0.1488], [0.1488, -0.1488]])
            qt = np.concatenate([qt, self.postprocessing([theta1, theta2, theta3])])
        return qt.transpose().tolist()


class IKSolverCUTER3DoFNum(_IKSolverCUTER3DoF):
    """
    Numerical IK
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        pass


class _IKSolverCUTER6DoF(_IKSolverCUTER):
    def __init__(self):
        super().__init__()
        self.theta_min = [-90, 0, -140, -180, -110, -180]
        self.theta_max = [90, 180, 45, 180, 110, 180]
        self.l1 = 10.18
        self.l2 = 19.41
        self.l3 = 2.91
        self.l4 = 25.22
        self.l5 = 3.00

    def postprocessing(self, rad):
        pass

    def __call__(self, taskspace: list):
        raise NotImplementedError


class IKSolverCUTER6DoFAna(_IKSolverCUTER6DoF):
    """
    Analytical IK
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        pass


class IKSolverCUTER6DoFNum(_IKSolverCUTER6DoF):
    """
    Numerical IK
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        pass
