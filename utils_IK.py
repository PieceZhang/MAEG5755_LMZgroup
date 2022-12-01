# Copyright (c) Zhang Yuelin. All Rights Reserved.
import numpy as np
from functools import partial
from sympy import Symbol, solve
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

    def __call__(self, taskspace: list):
        raise NotImplementedError


class IKSolverCUTER3DoFAna(_IKSolverCUTER):
    """
    Analytical IK
    """

    def __init__(self):
        super().__init__()

    def __call__(self, taskspace: list):
        pass


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
