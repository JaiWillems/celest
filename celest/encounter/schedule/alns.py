

from typing import Any, List, Callable
import math
import random


_BEST_SCORE = 3
_BETTER_SCORE = 2
_ACCEPT_SCORE = 1
_REJECT_SCORE = 0

Q = 1


class ALNS:
    """ALNS(x_init)

    The ANLS class is a generic framework of the adaptive large
    neighborhood search algorithm.

    Parameters
    ----------
    x_init : Any
        Initial solution.
    """

    def __init__(self, x_init: Any) -> None:

        self.x_init = x_init

    def add_cost_func(self, c: Callable) -> None:
        """Add cost function to minimize.

        Parameters
        ----------
        c : function
            Cost function that takes a solution as an input and returns a
            scalar cost value.
        """

        self.c = c

    def add_destroy_funcs(self, funcs: List[Callable]) -> None:
        """Add destroy functions to the ALNS instance.

        Parameters
        ----------
        funcs : List[function]
            List of destroy functions that take a solution as an input and
            return a destroyed solution.
        """

        self.destroy_funcs = funcs

    def add_repair_funcs(self, funcs: List[Callable]) -> None:
        """Add repair functions to the ALNS instance.

        Parameters
        ----------
        funcs : List[function]
            List of repair functions that take a destroyed solution as an input
            and return a repaired solution.
        """

        self.repair_funcs = funcs

    def _get_probabilities(self, weights: list) -> list:
        """Get probabilities from weights.

        Parameters
        ----------
        weights : list
            List of weights.

        Returns
        -------
        list
            list of probabilities.
        """

        s = sum(weights)
        p = [w / s for w in weights]

        return p

    def _get_destroy_index(self) -> int:
        """Return index of the destroy function to be used.

        Returns
        -------
        int
            Index of the destroy function to be used.
        """

        p = self._get_probabilities(self.destroy_weights)
        i = random.choices([j for j in range(len(self.destroy_funcs))], p)[0]

        return i

    def _get_repair_index(self) -> int:
        """Return index of the repair function to be used.

        Returns
        -------
        int
            Index of the repair function to be used.
        """

        p = self._get_probabilities(self.repair_weights)
        i = random.choices([j for j in range(len(self.repair_funcs))], p)[0]

        return i

    def _get_destroy_func(self, i) -> Callable:
        """Return destroy function at index i.

        Parameters
        ----------
        i : int
            Destroy function index.

        Returns
        -------
        function
        """

        return self.destroy_funcs[i]

    def _get_repair_func(self, i) -> Callable:
        """Return repair function at index i.

        Parameters
        ----------
        i : int
            Repair function index.

        Returns
        -------
        function
        """

        return self.repair_funcs[i]

    def _update_destroy_weights(self, l: float, score: int, i: int) -> None:
        """Update the destroy weights.

        Parameters
        ----------
        l : float
            Decay parameter within the range [0, 1].

            Determines how sensitive the weights are to changes from the
            performance of the destroy function.
        score : int
            Score update of the destroy function.
        i : int
            Destroy function index.
        """

        self.destroy_weights[i] = l * self.destroy_weights[i] + (1 - l) * score

    def _update_repair_weights(self, l: float, score: int, i: int) -> None:
        """Update the repair weights.

        Parameters
        ----------
        l : float
            Decay parameter within the range [0, 1].

            Determines how sensitive the weights are to changes from the
            performance of the repair function.
        score : int
            Score update of the repair function.
        i : int
            Repair function index.
        """

        self.repair_weights[i] = l * self.repair_weights[i] + (1 - l) * score

    def _accept(self, xt: Any, x: Any, T: float) -> bool:
        """Accpetance criterion.

        This method implements a simulated annealing acceptance criterion.

        Parameters
        ----------
        xt : Any
            Temporary solution.
        x : Any
            Current solution.
        T : float
            Temperature.

        Returns
        -------
        bool
            True if the solution is accepted, False otherwise.
        """

        if self.c(xt) <= self.c(x):
            return True
        else:
            p = math.exp(-(self.c(xt) - self.c(x)) / T)
            return random.choices([False, True], [1 - p, p])[0]
    
    def add_is_complete_func(self, func) -> None:
        """Function to check if the solution is complete.

        Parameters
        ----------
        func : Callable
            Function that takes a solution and returns a boolean.
        """

        self.is_complete = func

    def solve(self, max_iter: int, t0: float, p: float, l: float) -> Any:
        """Determine the optimal solution.

        Parameters
        ----------
        max_iter : int
            Maximum number of iterations.
        t0 : float
            Initial temperature.
        p : float
            Annealing coefficient.
        l : float
            Decay parameter within the range [0, 1].

        Returns
        -------
        Any
            Optimal solution.
        """

        x = self.x_init
        xb = self.x_init
        t = t0
        self.destroy_weights = [1] * len(self.destroy_funcs)
        self.repair_weights = [1] * len(self.repair_funcs)

        if self.is_complete(xb):
            return xb

        for _ in range(max_iter):

            i = self._get_destroy_index()
            j = self._get_repair_index()

            xt = self._get_repair_func(j)(self._get_destroy_func(i)(x, Q), Q)

            score = _REJECT_SCORE

            if self._accept(xt, x, t):
                x = xt
                score = _BETTER_SCORE if self.c(xt) > self.c(x) else _ACCEPT_SCORE

            if self.c(xt) < self.c(xb):
                xb = xt
                score = _BEST_SCORE

            if self.is_complete(xb):
                return xb

            self._update_destroy_weights(l, score, i)
            self._update_repair_weights(l, score, j)

            t = p * t

        return xb
