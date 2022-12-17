# Copyright (c) Zhang Yuelin. All Rights Reserved.
import math
import numpy as np
from functools import partial
from numpy import sin, cos, arccos, arcsin, arctan2, pi, sqrt
import matplotlib.pyplot as plt
from utils_FK import CUTER_FK_3DOF, CUTER_FK_6DOF, CUTER_FK_6DOFxyz
import utils_equations as eq


def rad2deg(x):
    return x / pi * 180


def deg2rad(x):
    return x / 180 * pi


class _IKSolverCUTER(object):
    """
    Basic IK class
    """

    def __init__(self):
        """
        Basic IK solver class
        """
        self.theta_min_deg = [-90, -15, -140, -180, -110, -180]
        self.theta_max_deg = [90, 180, 45, 180, 110, 180]
        self.l1 = 10.18
        self.l2 = 19.41
        self.l3 = 2.91
        self.l4 = 25.22
        self.l5 = 3.00

    def postprocessing_3dof(self, rad):
        """
        post processing, convert rad roots to angle
        :param rad: rad
        :return: angle
        """
        deg = list(map(lambda x: rad2deg(x), rad))
        constrain = [[[self.theta_min_deg[0] < deg[0][0, 0] < self.theta_max_deg[0],
                       self.theta_min_deg[0] < deg[0][0, 1] < self.theta_max_deg[0]],
                      [self.theta_min_deg[0] < deg[0][1, 0] < self.theta_max_deg[0],
                       self.theta_min_deg[0] < deg[0][1, 1] < self.theta_max_deg[0]]],
                     [[self.theta_min_deg[1] < deg[1][0, 0] < self.theta_max_deg[1],
                       self.theta_min_deg[1] < deg[1][0, 1] < self.theta_max_deg[1]],
                      [self.theta_min_deg[1] < deg[1][1, 0] < self.theta_max_deg[1],
                       self.theta_min_deg[1] < deg[1][1, 1] < self.theta_max_deg[1]]],
                     [[self.theta_min_deg[2] < deg[2][0] < self.theta_max_deg[2],
                       self.theta_min_deg[2] < deg[2][1] < self.theta_max_deg[2]]]]
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

    def solve3dofana(self, x, y, z):
        """
        3 dof analytical IK
        :param x: x coor
        :param y: y coor
        :param z: z coor
        :return: theta1, theta2, theta3
        """
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
        # postprocessing
        q = self.postprocessing_3dof([theta1, theta2, theta3])
        # FK (for verification)
        # x_fk, y_fk, z_fk = CUTER_FK_3DOF(self, deg2rad(np.array([q[0, 0], q[0, 1], q[0, 2]])))
        return q

    def postprocessing_6dof(self, rad):
        deg = list(map(lambda x: rad2deg(x), rad))
        deg[4] = 180 - deg[4]
        return np.array([[deg[0], deg[1], deg[2], deg[3], deg[4], deg[5]]])

    def solve6dofana(self, x, y, z, alpha, beta, gamma):
        """
        6 dof analytical IK
        :param x: x coor
        :param y: y coor
        :param z: z coor
        :param alpha: x rot
        :param beta: y rot
        :param gamma: z rot
        :return: theta1, theta2, theta3, theta4, theta5, theta6
        """
        # solving for theta123
        q123 = self.solve3dofana(x, y, z)
        theta1, theta2, theta3 = deg2rad(q123[0, 0]), deg2rad(q123[0, 1]), deg2rad(q123[0, 2])
        # R36
        R03 = eq.R03mat(theta1, theta2, theta3)
        R06 = eq.R06mat(theta1, theta2, theta3, deg2rad(alpha), deg2rad(beta), deg2rad(gamma))
        R36 = np.linalg.inv(R03) @ R06
        # solving for theta5
        theta5 = arccos(R36[1, 2])  # np.array([theta5, -theta5])
        sintheta5 = sin(theta5)
        if sintheta5 > 1e-8:
            # solving for theta4
            theta4 = arctan2(R36[2, 2] / sintheta5, R36[0, 2] / sintheta5)
            # solving for theta6
            theta6 = arctan2(-R36[1, 1] / sintheta5, R36[1, 0] / sintheta5)
        else:
            theta4 = theta6 = 0
        return self.postprocessing_6dof([theta1, theta2, theta3, theta4, theta5, theta6])

    def solvenum(self, taskspace, J, FK, qtlast, dof=3, zfunc=None):
        """
        numerical IK (3 dof / 6 dof)
        :param zfunc: (null space operation) function for z
        :param dof: DoF
        :param taskspace: list of x
        :param J: Jacobian func
        :param FK: FK func
        :param qtlast: initial q
        :return: qt (in array)
        """
        taskspace = np.array(taskspace).transpose()
        dxr = np.diff(taskspace.transpose()).transpose()  # dx reference
        dxr = np.concatenate([dxr, dxr[-1, :][None, :]])
        dt = 0.02  # time step 50Hz
        qt = np.ndarray((0, dof))  # qt
        et = np.ndarray((0, taskspace.shape[1]))  # et
        if taskspace.shape[1] == 3:
            Kp = np.array([25 for _ in range(3)])
        else:
            Kp = np.array([25 for _ in range(3)] + [0.05 for _ in range(3)])

        def offset(_):
            x = _.copy()
            x[1] = x[1] + 0.1488
            x[2] = x[2] - 0.1488
            return x

        for i, xyz in enumerate(taskspace):
            if taskspace.shape[1] == 3:  # [x y z]
                xyz = np.array([-xyz[0], -xyz[2], xyz[1]])
            else:  # [x y z alpha beta gamma]
                xyz = np.array([-xyz[0], -xyz[2], xyz[1], deg2rad(xyz[3]), deg2rad(xyz[4]), deg2rad(xyz[5])])
            dxc = dxr[i, :] + Kp * (xyz - FK(q=offset(qtlast[0])))  # dx control
            Jq = J(qtlast)
            invJq = np.linalg.pinv(Jq)
            if zfunc is None:
                qtlast = qtlast + dt * (invJq @ dxc)
            else:  # if use null space
                z = zfunc(qtlast)  # null space vector
                qtlast = qtlast + (dt * (invJq @ dxc)[:, None] + (np.eye(dof) - invJq @ Jq) @ z).transpose()
            # limit
            for j in range(qtlast.shape[1]):
                qtlast[0, j] = max(deg2rad(self.theta_min_deg[j]), min(deg2rad(self.theta_max_deg[j]), qtlast[0, j]))
            # cat
            qt = np.concatenate([qt, list(map(lambda x: rad2deg(offset(x)), qtlast))[0][None, :]])
            # for debug
            et = np.concatenate([et, (xyz - FK(q=offset(qtlast[0])))[None, :]])
            print(xyz, FK(q=deg2rad(qt[-1])))

        return qt


class _IKSolverCUTER3DoF(_IKSolverCUTER):
    def __init__(self):
        super().__init__()

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
            q = self.solve3dofana(x, y, z)
            qt = np.concatenate([qt, q])
        return qt.transpose().tolist()


class IKSolverCUTER3DoFNum(_IKSolverCUTER3DoF):
    """
    Numerical IK
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):

        def J(_):
            # redefine variables name
            l1 = self.l1
            l2 = sqrt(self.l2 ** 2 + self.l3 ** 2)
            l3 = self.l4
            theta1 = _[0, 0]
            theta2 = _[0, 1]
            theta3 = _[0, 2]
            return eq.J_3dof(theta1, theta2, theta3, l1, l2, l3)

        initq = self.solve3dofana(-taskspace[0][0], -taskspace[2][0], taskspace[1][0])
        initq = deg2rad(initq)
        qt = self.solvenum(taskspace, J, partial(CUTER_FK_3DOF, ik=self), initq)
        return qt.transpose().tolist()


class _IKSolverCUTER6DoF(_IKSolverCUTER):
    def __init__(self):
        super().__init__()
        self.l4 += self.l5

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
            q = self.solve6dofana(x, y, z, alpha, beta, gamma)
            qt = np.concatenate([qt, q])
        return qt.transpose().tolist()


class IKSolverCUTER6DoFNum(_IKSolverCUTER6DoF):
    """
    Numerical IK
    x=[x y z alpha beta gamma]
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        def J(_):
            theta1 = _[0, 0]
            theta2 = _[0, 1]
            theta3 = _[0, 2]
            theta4 = _[0, 3]
            theta5 = _[0, 4]
            theta6 = _[0, 5]
            return eq.J_6dof(theta1, theta2, theta3, theta4, theta5, theta6,
                             self.l1, self.l2, self.l3, self.l4, self.l5)

        initq = self.solve3dofana(-taskspace[0][0], -taskspace[2][0], taskspace[1][0])
        initq = deg2rad(initq)
        initq = np.concatenate([initq, np.array([[0, 0, 0]])], axis=1)

        qt = self.solvenum(taskspace, J, partial(CUTER_FK_6DOF, ik=self), initq, dof=6)
        return qt.transpose().tolist()


class IKSolverCUTER6DoFNumxyz(_IKSolverCUTER6DoF):
    """
    Numerical IK
    x=[x y z]
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        def J(_):
            theta1 = _[0, 0]
            theta2 = _[0, 1]
            theta3 = _[0, 2]
            theta4 = _[0, 3]
            theta5 = _[0, 4]
            theta6 = _[0, 5]
            return eq.J_6dof_xyz(theta1, theta2, theta3, theta4, theta5, theta6,
                                 self.l1, self.l2, self.l3, self.l4, self.l5)

        initq = self.solve3dofana(-taskspace[0][0], -taskspace[2][0], taskspace[1][0])
        initq = deg2rad(initq)
        initq = np.concatenate([initq, np.array([[0, 0, 0]])], axis=1)

        def zfunc(qtlast):
            qd = deg2rad(np.array([[np.mean([self.theta_max_deg[0], self.theta_min_deg[0]])],
                                   [np.mean([self.theta_max_deg[1], self.theta_min_deg[1]])],
                                   [np.mean([self.theta_max_deg[2], self.theta_min_deg[2]])],
                                   [np.mean([self.theta_max_deg[3], self.theta_min_deg[3]])],
                                   [np.mean([self.theta_max_deg[4], self.theta_min_deg[4]])],
                                   [np.mean([self.theta_max_deg[5], self.theta_min_deg[5]])]]))
            a = 0.02
            z = -a * (qtlast.transpose() - qd)
            return z

        qt = self.solvenum(taskspace, J, partial(CUTER_FK_6DOFxyz, ik=self), initq, dof=6, zfunc=zfunc)
        # qt = self.solvenum(taskspace, J, partial(CUTER_FK_6DOFxyz, ik=self), initq, dof=6, zfunc=None)
        return qt.transpose().tolist()
