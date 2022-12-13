# Copyright (c) Zhang Yuelin. All Rights Reserved.
import numpy as np
from numpy import sin, cos, sqrt


def CUTER_FK_3DOF(ik, q):
    theta1 = q[0]
    theta2 = q[1]
    theta3 = q[2]
    l1 = ik.l1
    l2 = sqrt(ik.l2 ** 2 + ik.l3 ** 2)
    l3 = ik.l4
    x = -sin(theta1) * (l2 * cos(theta2 - 0.1488) + l3 * cos(theta2 + theta3))
    y = cos(theta1) * (l2 * cos(theta2 - 0.1488) + l3 * cos(theta2 + theta3))
    z = l1 + l2 * sin(theta2 - 0.1488) + l3 * sin(theta2 + theta3)
    return x, y, z


def CUTER_FK_6DOF(ik, q):
    pass
