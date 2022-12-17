import numpy as np
from numpy import sin, cos


def R03mat(theta1, theta2, theta3):
    return np.array([[-cos(theta2 + theta3) * sin(theta1), sin(theta2 + theta3) * sin(theta1), cos(theta1)],
                     [cos(theta2 + theta3) * cos(theta1), -sin(theta2 + theta3) * cos(theta1), sin(theta1)],
                     [sin(theta2 + theta3), cos(theta2 + theta3), 0]])


def R06mat(theta1, theta2, theta3, theta4, theta5, theta6):
    return np.array([[cos(theta6) * (
            cos(theta5) * (cos(theta1) * sin(theta4) - cos(theta2 + theta3) * cos(theta4) * sin(theta1)) + sin(
        theta2 + theta3) * sin(theta1) * sin(theta5)) + sin(theta6) * (
                              cos(theta1) * cos(theta4) + cos(theta2 + theta3) * sin(theta1) * sin(theta4)),
                      cos(theta6) * (
                              cos(theta1) * cos(theta4) + cos(theta2 + theta3) * sin(theta1) * sin(theta4)) - sin(
                          theta6) * (cos(theta5) * (
                              cos(theta1) * sin(theta4) - cos(theta2 + theta3) * cos(theta4) * sin(theta1)) + sin(
                          theta2 + theta3) * sin(theta1) * sin(theta5)), sin(theta5) * (
                              cos(theta1) * sin(theta4) - cos(theta2 + theta3) * cos(theta4) * sin(theta1)) - sin(
            theta2 + theta3) * cos(theta5) * sin(theta1)],
                     [sin(theta6) * (
                             cos(theta4) * sin(theta1) - cos(theta2 + theta3) * cos(theta1) * sin(theta4)) + cos(
                         theta6) * (cos(theta5) * (
                             sin(theta1) * sin(theta4) + cos(theta2 + theta3) * cos(theta1) * cos(theta4)) - sin(
                         theta2 + theta3) * cos(theta1) * sin(theta5)), cos(theta6) * (
                              cos(theta4) * sin(theta1) - cos(theta2 + theta3) * cos(theta1) * sin(theta4)) - sin(
                         theta6) * (cos(theta5) * (
                             sin(theta1) * sin(theta4) + cos(theta2 + theta3) * cos(theta1) * cos(theta4)) - sin(
                         theta2 + theta3) * cos(theta1) * sin(theta5)), sin(theta5) * (
                              sin(theta1) * sin(theta4) + cos(theta2 + theta3) * cos(theta1) * cos(theta4)) + sin(
                         theta2 + theta3) * cos(theta1) * cos(theta5)],
                     [cos(theta6) * (cos(theta2 + theta3) * sin(theta5) + sin(theta2 + theta3) * cos(theta4) * cos(
                         theta5)) - sin(theta2 + theta3) * sin(theta4) * sin(theta6), - sin(theta6) * (
                              cos(theta2 + theta3) * sin(theta5) + sin(theta2 + theta3) * cos(theta4) * cos(
                          theta5)) - sin(theta2 + theta3) * cos(theta6) * sin(theta4),
                      sin(theta2 + theta3) * cos(theta4) * sin(theta5) - cos(theta2 + theta3) * cos(theta5)]])


def J_3dof(theta1, theta2, theta3, l1, l2, l3):
    return np.array([[-cos(theta1) * (l3 * cos(theta2 + theta3) + l2 * cos(0.1488 - theta2)),
                      sin(theta1) * (l3 * sin(theta2 + theta3) - l2 * sin(0.1488 - theta2)),
                      l3 * sin(theta2 + theta3) * sin(theta1)],
                     [-sin(theta1) * (l3 * cos(theta2 + theta3) + l2 * cos(0.1488 - theta2)),
                      -cos(theta1) * (l3 * sin(theta2 + theta3) - l2 * sin(0.1488 - theta2)),
                      -l3 * sin(theta2 + theta3) * cos(theta1)],
                     [0,
                      l3 * cos(theta2 + theta3) + l2 * cos(0.1488 - theta2),
                      l3 * cos(theta2 + theta3)]])


def J_6dof(theta1, theta2, theta3, theta4, theta5, theta6, l1, l2, l3, l4, l5):
    return np.array(
        [[l5 * cos(theta1) * cos(theta2) * sin(theta3) * sin(theta5) - cos(theta1) * cos(0.1488 - theta2) * (
                l2 ** 2 + l3 ** 2) ** (1 / 2) - l5 * cos(theta5) * sin(theta1) * sin(theta4) - l4 * cos(
            theta2 + theta3) * cos(theta1) + l5 * cos(theta1) * cos(theta3) * sin(theta2) * sin(theta5) - l5 * cos(
            theta1) * cos(theta2) * cos(theta3) * cos(theta4) * cos(theta5) + l5 * cos(theta1) * cos(theta4) * cos(
            theta5) * sin(theta2) * sin(theta3), sin(theta1) * (
                  l4 * sin(theta2 + theta3) - sin(0.1488 - theta2) * (l2 ** 2 + l3 ** 2) ** (
                  1 / 2) + l5 * cos(theta2) * cos(theta3) * sin(theta5) - l5 * sin(theta2) * sin(
              theta3) * sin(theta5) + l5 * cos(theta2) * cos(theta4) * cos(theta5) * sin(
              theta3) + l5 * cos(theta3) * cos(theta4) * cos(theta5) * sin(theta2)), sin(theta1) * (
                  l4 * sin(theta2 + theta3) + l5 * cos(theta2) * cos(theta3) * sin(theta5) - l5 * sin(
              theta2) * sin(theta3) * sin(theta5) + l5 * cos(theta2) * cos(theta4) * cos(theta5) * sin(
              theta3) + l5 * cos(theta3) * cos(theta4) * cos(theta5) * sin(theta2)),
          l5 * cos(theta5) * (cos(theta1) * cos(theta4) + cos(theta2) * cos(theta3) * sin(theta1) * sin(
              theta4) - sin(theta1) * sin(theta2) * sin(theta3) * sin(theta4)),
          l5 * cos(theta2) * cos(theta5) * sin(theta1) * sin(theta3) - l5 * cos(theta1) * sin(theta4) * sin(
              theta5) + l5 * cos(theta3) * cos(theta5) * sin(theta1) * sin(theta2) + l5 * cos(theta2) * cos(
              theta3) * cos(theta4) * sin(theta1) * sin(theta5) - l5 * cos(theta4) * sin(theta1) * sin(
              theta2) * sin(theta3) * sin(theta5), 0],
         [l5 * cos(theta1) * cos(theta5) * sin(theta4) - sin(theta1) * cos(0.1488 - theta2) * (
                 l2 ** 2 + l3 ** 2) ** (1 / 2) - l4 * cos(theta2 + theta3) * sin(theta1) + l5 * cos(
             theta2) * sin(theta1) * sin(theta3) * sin(theta5) + l5 * cos(theta3) * sin(theta1) * sin(
             theta2) * sin(theta5) - l5 * cos(theta2) * cos(theta3) * cos(theta4) * cos(theta5) * sin(
             theta1) + l5 * cos(theta4) * cos(theta5) * sin(theta1) * sin(theta2) * sin(theta3),
          -cos(theta1) * (l4 * sin(theta2 + theta3) - sin(0.1488 - theta2) * (l2 ** 2 + l3 ** 2) ** (
                  1 / 2) + l5 * cos(theta2) * cos(theta3) * sin(theta5) - l5 * sin(theta2) * sin(
              theta3) * sin(theta5) + l5 * cos(theta2) * cos(theta4) * cos(theta5) * sin(theta3) + l5 * cos(
              theta3) * cos(theta4) * cos(theta5) * sin(theta2)), -cos(theta1) * (
                  l4 * sin(theta2 + theta3) + l5 * cos(theta2) * cos(theta3) * sin(theta5) - l5 * sin(
              theta2) * sin(theta3) * sin(theta5) + l5 * cos(theta2) * cos(theta4) * cos(theta5) * sin(
              theta3) + l5 * cos(theta3) * cos(theta4) * cos(theta5) * sin(theta2)),
          l5 * cos(theta5) * (cos(theta4) * sin(theta1) - cos(theta1) * cos(theta2) * cos(theta3) * sin(
              theta4) + cos(theta1) * sin(theta2) * sin(theta3) * sin(theta4)),
          l5 * cos(theta1) * cos(theta4) * sin(theta2) * sin(theta3) * sin(theta5) - l5 * cos(theta1) * cos(
              theta2) * cos(theta5) * sin(theta3) - l5 * cos(theta1) * cos(theta3) * cos(theta5) * sin(
              theta2) - l5 * cos(theta1) * cos(theta2) * cos(theta3) * cos(theta4) * sin(theta5) - l5 * sin(
              theta1) * sin(theta4) * sin(theta5), 0],
         [0, l4 * cos(theta2 + theta3) + cos(0.1488 - theta2) * (l2 ** 2 + l3 ** 2) ** (1 / 2) - l5 * cos(
             theta2) * sin(theta3) * sin(theta5) - l5 * cos(theta3) * sin(theta2) * sin(theta5) + l5 * cos(
             theta2) * cos(theta3) * cos(theta4) * cos(theta5) - l5 * cos(theta4) * cos(theta5) * sin(
             theta2) * sin(theta3),
          l4 * cos(theta2) * cos(theta3) - l4 * sin(theta2) * sin(theta3) - l5 * cos(theta2) * sin(
              theta3) * sin(theta5) - l5 * cos(theta3) * sin(theta2) * sin(theta5) + l5 * cos(theta2) * cos(
              theta3) * cos(theta4) * cos(theta5) - l5 * cos(theta4) * cos(theta5) * sin(theta2) * sin(
              theta3), -l5 * sin(theta2 + theta3) * cos(theta5) * sin(theta4), -l5 * (
                  cos(theta5) * sin(theta2) * sin(theta3) - cos(theta2) * cos(theta3) * cos(
              theta5) + cos(theta2) * cos(theta4) * sin(theta3) * sin(theta5) + cos(theta3) * cos(
              theta4) * sin(theta2) * sin(theta5)), 0],
         [0, cos(theta1), cos(theta1), -sin(theta2 + theta3) * sin(theta1),
          cos(theta1) * cos(theta4) + cos(theta2) * cos(theta3) * sin(theta1) * sin(theta4) - sin(
              theta1) * sin(theta2) * sin(theta3) * sin(theta4),
          cos(theta1) * sin(theta4) * sin(theta5) - cos(theta2) * cos(theta5) * sin(theta1) * sin(
              theta3) - cos(theta3) * cos(theta5) * sin(theta1) * sin(theta2) - cos(theta2) * cos(
              theta3) * cos(theta4) * sin(theta1) * sin(theta5) + cos(theta4) * sin(theta1) * sin(
              theta2) * sin(theta3) * sin(theta5)],
         [0, sin(theta1), sin(theta1), sin(theta2 + theta3) * cos(theta1),
          cos(theta4) * sin(theta1) - cos(theta1) * cos(theta2) * cos(theta3) * sin(theta4) + cos(
              theta1) * sin(theta2) * sin(theta3) * sin(theta4),
          sin(theta1) * sin(theta4) * sin(theta5) + cos(theta1) * cos(theta2) * cos(theta5) * sin(
              theta3) + cos(theta1) * cos(theta3) * cos(theta5) * sin(theta2) + cos(theta1) * cos(
              theta2) * cos(theta3) * cos(theta4) * sin(theta5) - cos(theta1) * cos(theta4) * sin(
              theta2) * sin(theta3) * sin(theta5)],
         [1, 0, 0, -cos(theta2 + theta3), -sin(theta2 + theta3) * sin(theta4),
          cos(theta5) * sin(theta2) * sin(theta3) - cos(theta2) * cos(theta3) * cos(theta5) + cos(
              theta2) * cos(theta4) * sin(theta3) * sin(theta5) + cos(theta3) * cos(theta4) * sin(
              theta2) * sin(theta5)]])


# def J_6dof(theta1, theta2, theta3, theta4, theta5, theta6, l1, l2, l3, l4, l5):
#     return np.array([[l5 * cos(theta1) * cos(theta4) * sin(theta2) * sin(theta3) * sin(theta5) - sin(0.1488) * cos(
#         theta1) * sin(theta2) * (l2 ** 2 + l3 ** 2) ** (1 / 2) - l4 * cos(theta1) * cos(theta2) * sin(theta3) - l4 * cos(
#         theta1) * cos(theta3) * sin(theta2) - l5 * sin(theta1) * sin(theta4) * sin(theta5) - l5 * cos(theta1) * cos(
#         theta2) * cos(theta5) * sin(theta3) - l5 * cos(theta1) * cos(theta3) * cos(theta5) * sin(theta2) - l5 * cos(
#         theta1) * cos(theta2) * cos(theta3) * cos(theta4) * sin(theta5) - cos(0.1488) * cos(theta1) * cos(theta2) * (
#                                   l2 ** 2 + l3 ** 2) ** (1 / 2), sin(theta1) * (
#                                   cos(0.1488) * sin(theta2) * (l2 ** 2 + l3 ** 2) ** (1 / 2) - sin(0.1488) * cos(theta2) * (
#                                       l2 ** 2 + l3 ** 2) ** (1 / 2) - l4 * cos(theta2) * cos(theta3) + l4 * sin(
#                               theta2) * sin(theta3) - l5 * cos(theta2) * cos(theta3) * cos(theta5) + l5 * cos(
#                               theta5) * sin(theta2) * sin(theta3) + l5 * cos(theta2) * cos(theta4) * sin(theta3) * sin(
#                               theta5) + l5 * cos(theta3) * cos(theta4) * sin(theta2) * sin(theta5)), sin(theta1) * (
#                                   cos(0.1488) * sin(theta2) * (l2 ** 2 + l3 ** 2) ** (1 / 2) - sin(0.1488) * cos(theta2) * (
#                                       l2 ** 2 + l3 ** 2) ** (1 / 2) - l4 * cos(theta2) * cos(theta3) + l4 * sin(
#                               theta2) * sin(theta3) - l5 * cos(theta2) * cos(theta3) * cos(theta5) + l5 * cos(
#                               theta5) * sin(theta2) * sin(theta3) + l5 * cos(theta2) * cos(theta4) * sin(theta3) * sin(
#                               theta5) + l5 * cos(theta3) * cos(theta4) * sin(theta2) * sin(theta5)),
#                       l5 * sin(theta5) * (cos(theta1) * cos(theta4) + cos(theta2) * cos(theta3) * sin(theta1) * sin(
#                           theta4) - sin(theta1) * sin(theta2) * sin(theta3) * sin(theta4)),
#                       l5 * cos(theta1) * cos(theta5) * sin(theta4) + l5 * cos(theta2) * sin(theta1) * sin(theta3) * sin(
#                           theta5) + l5 * cos(theta3) * sin(theta1) * sin(theta2) * sin(theta5) - l5 * cos(theta2) * cos(
#                           theta3) * cos(theta4) * cos(theta5) * sin(theta1) + l5 * cos(theta4) * cos(theta5) * sin(
#                           theta1) * sin(theta2) * sin(theta3), 0],
#                      [l5 * cos(theta1) * sin(theta4) * sin(theta5) - sin(0.1488) * sin(theta1) * sin(theta2) * (
#                                  l2 ** 2 + l3 ** 2) ** (1 / 2) - l4 * cos(theta2) * sin(theta1) * sin(theta3) - l4 * cos(
#                          theta3) * sin(theta1) * sin(theta2) - cos(0.1488) * cos(theta2) * sin(theta1) * (
#                                   l2 ** 2 + l3 ** 2) ** (1 / 2) - l5 * cos(theta2) * cos(theta5) * sin(theta1) * sin(
#                          theta3) - l5 * cos(theta3) * cos(theta5) * sin(theta1) * sin(theta2) - l5 * cos(theta2) * cos(
#                          theta3) * cos(theta4) * sin(theta1) * sin(theta5) + l5 * cos(theta4) * sin(theta1) * sin(
#                          theta2) * sin(theta3) * sin(theta5), -cos(theta1) * (
#                                   cos(0.1488) * sin(theta2) * (l2 ** 2 + l3 ** 2) ** (1 / 2) - sin(0.1488) * cos(theta2) * (
#                                       l2 ** 2 + l3 ** 2) ** (1 / 2) - l4 * cos(theta2) * cos(theta3) + l4 * sin(
#                               theta2) * sin(theta3) - l5 * cos(theta2) * cos(theta3) * cos(theta5) + l5 * cos(
#                               theta5) * sin(theta2) * sin(theta3) + l5 * cos(theta2) * cos(theta4) * sin(theta3) * sin(
#                               theta5) + l5 * cos(theta3) * cos(theta4) * sin(theta2) * sin(theta5)), -cos(theta1) * (
#                                   cos(0.1488) * sin(theta2) * (l2 ** 2 + l3 ** 2) ** (1 / 2) - sin(0.1488) * cos(theta2) * (
#                                       l2 ** 2 + l3 ** 2) ** (1 / 2) - l4 * cos(theta2) * cos(theta3) + l4 * sin(
#                               theta2) * sin(theta3) - l5 * cos(theta2) * cos(theta3) * cos(theta5) + l5 * cos(
#                               theta5) * sin(theta2) * sin(theta3) + l5 * cos(theta2) * cos(theta4) * sin(theta3) * sin(
#                               theta5) + l5 * cos(theta3) * cos(theta4) * sin(theta2) * sin(theta5)),
#                       l5 * sin(theta5) * (cos(theta4) * sin(theta1) - cos(theta1) * cos(theta2) * cos(theta3) * sin(
#                           theta4) + cos(theta1) * sin(theta2) * sin(theta3) * sin(theta4)),
#                       l5 * cos(theta5) * sin(theta1) * sin(theta4) - l5 * cos(theta1) * cos(theta2) * sin(theta3) * sin(
#                           theta5) - l5 * cos(theta1) * cos(theta3) * sin(theta2) * sin(theta5) + l5 * cos(theta1) * cos(
#                           theta2) * cos(theta3) * cos(theta4) * cos(theta5) - l5 * cos(theta1) * cos(theta4) * cos(
#                           theta5) * sin(theta2) * sin(theta3), 0],
#                      [0, cos(0.1488) * cos(theta2) * (l2 ** 2 + l3 ** 2) ** (1 / 2) + sin(0.1488) * sin(theta2) * (
#                                  l2 ** 2 + l3 ** 2) ** (1 / 2) + l4 * cos(theta2) * sin(theta3) + l4 * cos(theta3) * sin(
#                          theta2) + l5 * cos(theta2) * cos(theta5) * sin(theta3) + l5 * cos(theta3) * cos(theta5) * sin(
#                          theta2) + l5 * cos(theta2) * cos(theta3) * cos(theta4) * sin(theta5) - l5 * cos(theta4) * sin(
#                          theta2) * sin(theta3) * sin(theta5),
#                       cos(0.1488) * cos(theta2) * (l2 ** 2 + l3 ** 2) ** (1 / 2) + sin(0.1488) * sin(theta2) * (
#                                   l2 ** 2 + l3 ** 2) ** (1 / 2) + l4 * cos(theta2) * sin(theta3) + l4 * cos(theta3) * sin(
#                           theta2) + l5 * cos(theta2) * cos(theta5) * sin(theta3) + l5 * cos(theta3) * cos(theta5) * sin(
#                           theta2) + l5 * cos(theta2) * cos(theta3) * cos(theta4) * sin(theta5) - l5 * cos(theta4) * sin(
#                           theta2) * sin(theta3) * sin(theta5), -l5 * sin(theta2 + theta3) * sin(theta4) * sin(theta5),
#                       l5 * cos(theta2) * cos(theta3) * sin(theta5) - l5 * sin(theta2) * sin(theta3) * sin(
#                           theta5) + l5 * cos(theta2) * cos(theta4) * cos(theta5) * sin(theta3) + l5 * cos(theta3) * cos(
#                           theta4) * cos(theta5) * sin(theta2), 0],
#                      [0, cos(theta1), cos(theta1), -sin(theta2 + theta3) * sin(theta1),
#                       cos(theta1) * cos(theta4) + cos(theta2) * cos(theta3) * sin(theta1) * sin(theta4) - sin(
#                           theta1) * sin(theta2) * sin(theta3) * sin(theta4),
#                       cos(theta1) * sin(theta4) * sin(theta5) - cos(theta2) * cos(theta5) * sin(theta1) * sin(
#                           theta3) - cos(theta3) * cos(theta5) * sin(theta1) * sin(theta2) - cos(theta2) * cos(
#                           theta3) * cos(theta4) * sin(theta1) * sin(theta5) + cos(theta4) * sin(theta1) * sin(
#                           theta2) * sin(theta3) * sin(theta5)],
#                      [0, sin(theta1), sin(theta1), sin(theta2 + theta3) * cos(theta1),
#                       cos(theta4) * sin(theta1) - cos(theta1) * cos(theta2) * cos(theta3) * sin(theta4) + cos(
#                           theta1) * sin(theta2) * sin(theta3) * sin(theta4),
#                       sin(theta1) * sin(theta4) * sin(theta5) + cos(theta1) * cos(theta2) * cos(theta5) * sin(
#                           theta3) + cos(theta1) * cos(theta3) * cos(theta5) * sin(theta2) + cos(theta1) * cos(
#                           theta2) * cos(theta3) * cos(theta4) * sin(theta5) - cos(theta1) * cos(theta4) * sin(
#                           theta2) * sin(theta3) * sin(theta5)],
#                      [1, 0, 0, -cos(theta2 + theta3), -sin(theta2 + theta3) * sin(theta4),
#                       cos(theta5) * sin(theta2) * sin(theta3) - cos(theta2) * cos(theta3) * cos(theta5) + cos(
#                           theta2) * cos(theta4) * sin(theta3) * sin(theta5) + cos(theta3) * cos(theta4) * sin(
#                           theta2) * sin(theta5)]])


def J_6dof_xyz(theta1, theta2, theta3, theta4, theta5, theta6, l1, l2, l3, l4, l5):
    return np.array(
        [[l5 * cos(theta1) * cos(theta2) * sin(theta3) * sin(theta5) - cos(theta1) * cos(0.1488 - theta2) * (
                l2 ** 2 + l3 ** 2) ** (1 / 2) - l5 * cos(theta5) * sin(theta1) * sin(theta4) - l4 * cos(
            theta2 + theta3) * cos(theta1) + l5 * cos(theta1) * cos(theta3) * sin(theta2) * sin(theta5) - l5 * cos(
            theta1) * cos(theta2) * cos(theta3) * cos(theta4) * cos(theta5) + l5 * cos(theta1) * cos(theta4) * cos(
            theta5) * sin(theta2) * sin(theta3), sin(theta1) * (
                  l4 * sin(theta2 + theta3) - sin(0.1488 - theta2) * (l2 ** 2 + l3 ** 2) ** (
                  1 / 2) + l5 * cos(theta2) * cos(theta3) * sin(theta5) - l5 * sin(theta2) * sin(
              theta3) * sin(theta5) + l5 * cos(theta2) * cos(theta4) * cos(theta5) * sin(
              theta3) + l5 * cos(theta3) * cos(theta4) * cos(theta5) * sin(theta2)), sin(theta1) * (
                  l4 * sin(theta2 + theta3) + l5 * cos(theta2) * cos(theta3) * sin(theta5) - l5 * sin(
              theta2) * sin(theta3) * sin(theta5) + l5 * cos(theta2) * cos(theta4) * cos(theta5) * sin(
              theta3) + l5 * cos(theta3) * cos(theta4) * cos(theta5) * sin(theta2)),
          l5 * cos(theta5) * (cos(theta1) * cos(theta4) + cos(theta2) * cos(theta3) * sin(theta1) * sin(
              theta4) - sin(theta1) * sin(theta2) * sin(theta3) * sin(theta4)),
          l5 * cos(theta2) * cos(theta5) * sin(theta1) * sin(theta3) - l5 * cos(theta1) * sin(theta4) * sin(
              theta5) + l5 * cos(theta3) * cos(theta5) * sin(theta1) * sin(theta2) + l5 * cos(theta2) * cos(
              theta3) * cos(theta4) * sin(theta1) * sin(theta5) - l5 * cos(theta4) * sin(theta1) * sin(
              theta2) * sin(theta3) * sin(theta5), 0],
         [l5 * cos(theta1) * cos(theta5) * sin(theta4) - sin(theta1) * cos(0.1488 - theta2) * (
                 l2 ** 2 + l3 ** 2) ** (1 / 2) - l4 * cos(theta2 + theta3) * sin(theta1) + l5 * cos(
             theta2) * sin(theta1) * sin(theta3) * sin(theta5) + l5 * cos(theta3) * sin(theta1) * sin(
             theta2) * sin(theta5) - l5 * cos(theta2) * cos(theta3) * cos(theta4) * cos(theta5) * sin(
             theta1) + l5 * cos(theta4) * cos(theta5) * sin(theta1) * sin(theta2) * sin(theta3),
          -cos(theta1) * (l4 * sin(theta2 + theta3) - sin(0.1488 - theta2) * (l2 ** 2 + l3 ** 2) ** (
                  1 / 2) + l5 * cos(theta2) * cos(theta3) * sin(theta5) - l5 * sin(theta2) * sin(
              theta3) * sin(theta5) + l5 * cos(theta2) * cos(theta4) * cos(theta5) * sin(theta3) + l5 * cos(
              theta3) * cos(theta4) * cos(theta5) * sin(theta2)), -cos(theta1) * (
                  l4 * sin(theta2 + theta3) + l5 * cos(theta2) * cos(theta3) * sin(theta5) - l5 * sin(
              theta2) * sin(theta3) * sin(theta5) + l5 * cos(theta2) * cos(theta4) * cos(theta5) * sin(
              theta3) + l5 * cos(theta3) * cos(theta4) * cos(theta5) * sin(theta2)),
          l5 * cos(theta5) * (cos(theta4) * sin(theta1) - cos(theta1) * cos(theta2) * cos(theta3) * sin(
              theta4) + cos(theta1) * sin(theta2) * sin(theta3) * sin(theta4)),
          l5 * cos(theta1) * cos(theta4) * sin(theta2) * sin(theta3) * sin(theta5) - l5 * cos(theta1) * cos(
              theta2) * cos(theta5) * sin(theta3) - l5 * cos(theta1) * cos(theta3) * cos(theta5) * sin(
              theta2) - l5 * cos(theta1) * cos(theta2) * cos(theta3) * cos(theta4) * sin(theta5) - l5 * sin(
              theta1) * sin(theta4) * sin(theta5), 0],
         [0, l4 * cos(theta2 + theta3) + cos(0.1488 - theta2) * (l2 ** 2 + l3 ** 2) ** (1 / 2) - l5 * cos(
             theta2) * sin(theta3) * sin(theta5) - l5 * cos(theta3) * sin(theta2) * sin(theta5) + l5 * cos(
             theta2) * cos(theta3) * cos(theta4) * cos(theta5) - l5 * cos(theta4) * cos(theta5) * sin(
             theta2) * sin(theta3),
          l4 * cos(theta2) * cos(theta3) - l4 * sin(theta2) * sin(theta3) - l5 * cos(theta2) * sin(
              theta3) * sin(theta5) - l5 * cos(theta3) * sin(theta2) * sin(theta5) + l5 * cos(theta2) * cos(
              theta3) * cos(theta4) * cos(theta5) - l5 * cos(theta4) * cos(theta5) * sin(theta2) * sin(
              theta3), -l5 * sin(theta2 + theta3) * cos(theta5) * sin(theta4), -l5 * (
                  cos(theta5) * sin(theta2) * sin(theta3) - cos(theta2) * cos(theta3) * cos(
              theta5) + cos(theta2) * cos(theta4) * sin(theta3) * sin(theta5) + cos(theta3) * cos(
              theta4) * sin(theta2) * sin(theta5)), 0]])
