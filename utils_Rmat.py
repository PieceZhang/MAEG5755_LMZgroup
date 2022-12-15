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
