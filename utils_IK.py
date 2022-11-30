# Copyright (c) Zhang Yuelin. All Rights Reserved.
import numpy as np
from functools import partial
from sympy import Symbol, solve
import matplotlib.pyplot as plt


class IKSolver(object):
    def __init__(self, num_dofs=3):
        """
        Basic IK solver class
        :param num_dofs: number of DoFs (q, task space)
        """
        self.num_dofs = num_dofs

    def __call__(self, taskspace: list):
        raise NotImplementedError


class IKSolverAnalytical(IKSolver):
    def __init__(self, num_dofs=3):
        super().__init__(num_dofs=num_dofs)

    def __call__(self, taskspace: list):
        pass


class IKSolverNumerical(IKSolver):
    def __init__(self, num_dofs=3):
        super().__init__(num_dofs=num_dofs)

    def __call__(self, taskspace: list):
        pass
