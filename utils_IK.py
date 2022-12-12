# Copyright (c) Zhang Yuelin. All Rights Reserved.
import math
import numpy as np
from functools import partial
from numpy import sin, cos, arccos, arcsin, arctan2, pi, sqrt
# from sympy import Symbol, solve, nsolve, sin, cos, acos, atan, pi
import matplotlib.pyplot as plt
from utils_FK import CUTER_FK_3DOF, CUTER_FK_6DOF


def rad2deg(x):
    return x / pi * 180


def deg2rad(x):
    return x / 180 * pi


class _IKSolverCUTER(object):
    def __init__(self):
        """
        Basic IK solver class
        """
        self.theta_min = [-90, 0, -140, -180, -110, -180]
        self.theta_max = [90, 180, 45, 180, 110, 180]
        self.l1 = 10.18
        self.l2 = 19.41
        self.l3 = 2.91
        self.l4 = 25.22
        self.l5 = 3.00

    def solve3dofana(self, x, y, z):
        # redefine variables name
        l1 = self.l1
        l2 = sqrt(self.l2 ** 2 + self.l3 ** 2)
        l3 = self.l4
        # solving for theta3
        theta3 = arccos(max(-1, min(1, (x ** 2 + y ** 2 + (z - l1) ** 2 - l2 ** 2 - l3 ** 2) / (2 * l2 * l3))))
        theta3 = np.array([theta3, -theta3])
        # solving for theta2
        a = l2 + l3 * cos(theta3)
        b = l3 * sin(theta3)
        alpha = arctan2(b / sqrt(a ** 2 + b ** 2), a / sqrt(a ** 2 + b ** 2))
        asi = np.array([arcsin(max(-1, min(1, (z - l1) / sqrt(a[0] ** 2 + b[0] ** 2)))),
                        arcsin(max(-1, min(1, (z - l1) / sqrt(a[1] ** 2 + b[1] ** 2))))])
        theta2 = np.array([[asi[0] - alpha[0], pi - (asi[0] - alpha[0])],
                           [asi[1] - alpha[1], pi - (asi[1] - alpha[1])]])
        # solving for theta1
        sin1_0 = -x / (l2 * cos(theta2[0]) + l3 * cos(theta2[0] + theta3[0]))
        cos1_0 = y / (l2 * cos(theta2[0]) + l3 * cos(theta2[0] + theta3[0]))
        sin1_1 = -x / (l2 * cos(theta2[1]) + l3 * cos(theta2[1] + theta3[1]))
        cos1_1 = y / (l2 * cos(theta2[1]) + l3 * cos(theta2[1] + theta3[1]))
        theta1 = np.array([arctan2(sin1_0, cos1_0), arctan2(sin1_1, cos1_1)])
        # cat
        theta3 -= 0.1488
        theta2 += 0.1488
        # theta2 += np.array([[0.1488, -0.1488], [0.1488, -0.1488]])
        # FK (for verification)
        x_fk, y_fk, z_fk = CUTER_FK_3DOF(self, [theta1, theta2, theta3])
        return theta1, theta2, theta3

    def solvenum(self, taskspace, J, FK):
        taskspace = np.array(taskspace).transpose()
        dxr = np.diff(taskspace.transpose()).transpose()  # dx reference
        dxr = np.concatenate([dxr, dxr[-1, :][None, :]])
        dt = 1  # time step 50Hz ??
        qtlast = np.array([[0 for _ in range(taskspace.shape[1])]])
        qt = np.ndarray((0, taskspace.shape[1]))  # qt
        et = np.ndarray((0, taskspace.shape[1]))  # et
        Kp = np.array([0.05, 0.05, 0.05])

        def offset(_):
            x = _.copy()
            x[1] = x[1] + 0.1488
            x[2] = x[2] - 0.1488
            return x

        for i, xyz in enumerate(taskspace):
            xyz = np.array([-xyz[0], -xyz[2], xyz[1]])
            dxc = dxr[i, :] + Kp * (xyz - FK(q=offset(qtlast[0])))  # dx control
            qtlast = qtlast + dt * (np.linalg.inv(J(qtlast)) @ dxc)
            # cat
            qt = np.concatenate([qt, list(map(lambda x: rad2deg(offset(x)), qtlast))[0][None, :]])
            # limit
            qt[-1, 0] = max(self.theta_min[0], min(self.theta_max[0], qt[-1, 0]))
            qt[-1, 1] = max(self.theta_min[1], min(self.theta_max[1], qt[-1, 1]))
            qt[-1, 2] = max(self.theta_min[2], min(self.theta_max[2], qt[-1, 2]))
            # for debug
            et = np.concatenate([et, (xyz - FK(q=offset(qtlast[0])))[None, :]])
            print(xyz, FK(q=deg2rad(qt[-1])))

        return qt

    def solve6dofana(self, x, y, z, alpha, beta, gamma):
        theta1, theta2, theta3 = self.solve3dofana(x, y, z)
        R03 = np.array([[-cos(theta2 + theta3) * sin(theta1), sin(theta2 + theta3) * sin(theta1), cos(theta1)],
                        [cos(theta2 + theta3) * cos(theta1), -sin(theta2 + theta3) * cos(theta1), sin(theta1)],
                        [sin(theta2 + theta3), cos(theta2 + theta3), 0]])
        # TODO
        theta4 = None
        theta5 = None
        theta6 = None
        return theta1, theta2, theta3, theta4, theta5, theta6


class _IKSolverCUTER3DoF(_IKSolverCUTER):
    def __init__(self):
        super().__init__()
        # some variables in basic class are overrided
        self.theta_min = [-90, -15, -140]
        self.theta_max = [90, 180, 45]
        # self.l4 = 20.2
        self.l4 = 25.22  # if use grabber, then set l4 to be longer

    def postprocessing(self, rad):
        """
        post processing, convert rad roots to angle
        :param rad: rad
        :return: angle
        """
        deg = list(map(lambda x: rad2deg(x), rad))
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
        select = []
        for theta3 in range(2):
            # theta3
            if not constrain[2][0][theta3]:
                continue
            # theta2
            for theta2 in range(2):
                # theta1
                if not constrain[1][theta3][theta2] or not constrain[0][theta3][theta2]:
                    continue
                else:
                    theta1 = theta2
                    select.append([theta1, theta2, theta3])
                    break
        # if len(select) == 2:
        #     theta1, theta2, theta3 = select[0]  # to choose 1 or 0 if multi solutions exist
        if not constrain[2][0][theta3] or not constrain[1][theta3][theta2] or not ('theta1' in locals().keys()):
            raise ValueError('[IKSolver] Location cannot reach!')
        elif False in [constrain[2][0][theta3], constrain[1][theta3][theta2], constrain[0][theta3][theta1]]:
            raise ValueError('[IKSolver] Runtime error!')
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
        taskspace = np.array(taskspace).transpose()
        qt = np.ndarray((0, 3))
        for xyz in taskspace:
            x = -xyz[0]
            z = xyz[1]
            y = -xyz[2]  # exchange y and z to match the coor in the simulator
            theta1, theta2, theta3 = self.solve3dofana(x, y, z)
            qt = np.concatenate([qt, self.postprocessing([theta1, theta2, theta3])])
        return qt.transpose().tolist()


class IKSolverCUTER3DoFNum(_IKSolverCUTER3DoF):
    """
    Numerical IK
    """

    def __init__(self):
        super().__init__()

    def postprocessing(self, qt):
        # for t in range(qt.shape[0]):
        #     qt[t] = qt[t] % (np.sign(qt[t]) * 360)
        return qt

    def __call__(self, taskspace: list):
        # redefine variables name
        l1 = self.l1
        l2 = sqrt(self.l2 ** 2 + self.l3 ** 2)
        l3 = self.l4
        beta = 0.1488

        def J(_):
            theta1 = _[0, 0]
            theta2 = _[0, 1]
            theta3 = _[0, 2]
            return np.array([[-cos(theta1) * (l3 * cos(theta2 + theta3) + l2 * cos(beta - theta2)),
                              sin(theta1) * (l3 * sin(theta2 + theta3) - l2 * sin(beta - theta2)),
                              l3 * sin(theta2 + theta3) * sin(theta1)],
                             [-sin(theta1) * (l3 * cos(theta2 + theta3) + l2 * cos(beta - theta2)),
                              -cos(theta1) * (l3 * sin(theta2 + theta3) - l2 * sin(beta - theta2)),
                              -l3 * sin(theta2 + theta3) * cos(theta1)],
                             [0,
                              l3 * cos(theta2 + theta3) + l2 * cos(beta - theta2),
                              l3 * cos(theta2 + theta3)]])

        qt = self.solvenum(taskspace, J, partial(CUTER_FK_3DOF, ik=self))
        qt = self.postprocessing(qt)
        return qt.transpose().tolist()


class _IKSolverCUTER6DoF(_IKSolverCUTER):
    def __init__(self):
        super().__init__()
        self.theta_min = [-90, 0, -140, -180, -110, -180]
        self.theta_max = [90, 180, 45, 180, 110, 180]

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
        taskspace = np.array(taskspace).transpose()
        qt = np.ndarray((0, 6))
        for ts in taskspace:
            x = -ts[0]
            z = ts[1]
            y = -ts[2]  # exchange y and z to match the coor in the simulator
            alpha = ts[3]
            beta = ts[4]
            gamma = ts[5]
            theta1, theta2, theta3, theta4, theta5, theta6 = self.solve6dofana(x, y, z, alpha, beta, gamma)
            qt = np.concatenate([qt, self.postprocessing([theta1, theta2, theta3, theta4, theta5, theta6])])
        return qt.transpose().tolist()


class IKSolverCUTER6DoFNum(_IKSolverCUTER6DoF):
    """
    Numerical IK
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        pass
