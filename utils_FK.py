# Copyright (c) Zhang Yuelin. All Rights Reserved.
import numpy as np
from numpy import sin, cos, sqrt, arcsin, arctan2


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


def CUTER_FK_6DOFxyz(ik, q):
    theta1 = q[0]
    theta2 = q[1]
    theta3 = q[2]
    theta4 = q[3]
    theta5 = q[4]
    theta6 = q[5]
    l1 = ik.l1
    l2 = ik.l2
    l3 = ik.l3
    l4 = ik.l4
    l5 = ik.l5
    beta = 0.1488
    x = l5 * cos(theta1) * cos(theta5) * sin(theta4) - sin(theta1) * cos(beta - theta2) * (l2 ** 2 + l3 ** 2) ** (
            1 / 2) - l4 * cos(theta2 + theta3) * sin(theta1) + l5 * cos(theta2) * sin(theta1) * sin(theta3) * sin(
        theta5) + l5 * cos(theta3) * sin(theta1) * sin(theta2) * sin(theta5) - l5 * cos(theta2) * cos(theta3) * cos(
        theta4) * cos(theta5) * sin(theta1) + l5 * cos(theta4) * cos(theta5) * sin(theta1) * sin(theta2) * sin(theta3)

    y = l4 * cos(theta2 + theta3) * cos(theta1) + cos(theta1) * cos(beta - theta2) * (l2 ** 2 + l3 ** 2) ** (
            1 / 2) + l5 * cos(theta5) * sin(theta1) * sin(theta4) - l5 * cos(theta1) * cos(theta2) * sin(
        theta3) * sin(theta5) - l5 * cos(theta1) * cos(theta3) * sin(theta2) * sin(theta5) + l5 * cos(
        theta1) * cos(theta2) * cos(theta3) * cos(theta4) * cos(theta5) - l5 * cos(theta1) * cos(theta4) * cos(
        theta5) * sin(theta2) * sin(theta3)

    z = l1 - l4 * (sin(beta - theta2) * cos(beta + theta3) - sin(beta + theta3) * cos(beta - theta2)) + l5 * (
            sin(theta5) * (
            sin(beta - theta2) * sin(beta + theta3) + cos(beta + theta3) * cos(beta - theta2)) - cos(
        theta4) * cos(theta5) * (sin(beta - theta2) * cos(beta + theta3) - sin(beta + theta3) * cos(
        beta - theta2))) - sin(beta - theta2) * (l2 ** 2 + l3 ** 2) ** (1 / 2)

    return x, y, z


def CUTER_FK_6DOF(ik, q):
    theta1 = q[0]
    theta2 = q[1]
    theta3 = q[2]
    theta4 = q[3]
    theta5 = q[4]
    theta6 = q[5]

    x, y, z = CUTER_FK_6DOFxyz(ik, q)

    beta = arcsin(sin(theta5) * (cos(theta1) * sin(theta4) - cos(theta4) * (
            cos(0.1488 + theta3) * sin(theta1) * cos(0.1488 - theta2) + sin(0.1488 - theta2) * sin(
        0.1488 + theta3) * sin(
        theta1))) + cos(theta5) * (
                          sin(0.1488 - theta2) * cos(0.1488 + theta3) * sin(theta1) - sin(0.1488 + theta3) * sin(
                      theta1) * cos(0.1488 - theta2)))

    gamma = arctan2(-(cos(theta6) * (cos(theta1) * cos(theta4) + sin(theta4) * (
            cos(0.1488 + theta3) * sin(theta1) * cos(0.1488 - theta2) + sin(0.1488 - theta2) * sin(
        0.1488 + theta3) * sin(
        theta1))) - sin(theta6) * (cos(theta5) * (cos(theta1) * sin(theta4) - cos(theta4) * (
            cos(0.1488 + theta3) * sin(theta1) * cos(0.1488 - theta2) + sin(0.1488 - theta2) * sin(
        0.1488 + theta3) * sin(
        theta1))) - sin(theta5) * (sin(0.1488 - theta2) * cos(0.1488 + theta3) * sin(theta1) - sin(
        0.1488 + theta3) * sin(
        theta1) * cos(0.1488 - theta2)))) / cos(beta),
                    (sin(theta6) * (cos(theta1) * cos(theta4) + sin(theta4) * (
                            cos(0.1488 + theta3) * sin(theta1) * cos(0.1488 - theta2) + sin(0.1488 - theta2) * sin(
                        0.1488 + theta3) * sin(theta1))) + cos(theta6) * (cos(theta5) * (
                            cos(theta1) * sin(theta4) - cos(theta4) * (
                            cos(0.1488 + theta3) * sin(theta1) * cos(0.1488 - theta2) + sin(0.1488 - theta2) * sin(
                        0.1488 + theta3) * sin(theta1))) - sin(theta5) * (sin(0.1488 - theta2) * cos(
                        0.1488 + theta3) * sin(theta1) - sin(0.1488 + theta3) * sin(theta1) * cos(
                        0.1488 - theta2)))) / cos(
                        beta))

    alpha = arctan2(-(sin(theta5) * (sin(theta1) * sin(theta4) + cos(theta4) * (
            cos(0.1488 + theta3) * cos(theta1) * cos(0.1488 - theta2) + sin(0.1488 - theta2) * sin(
        0.1488 + theta3) * cos(
        theta1))) - cos(theta5) * (sin(0.1488 - theta2) * cos(0.1488 + theta3) * cos(theta1) - sin(
        0.1488 + theta3) * cos(
        theta1) * cos(0.1488 - theta2))) / cos(beta),
                    (- cos(theta5) * (sin(0.1488 - theta2) * sin(0.1488 + theta3) + cos(0.1488 + theta3) * cos(
                        0.1488 - theta2)) - cos(theta4) * sin(theta5) * (
                             sin(0.1488 - theta2) * cos(0.1488 + theta3) - sin(0.1488 + theta3) * cos(
                         0.1488 - theta2))) / cos(beta))

    return x, y, z, alpha, beta, gamma
