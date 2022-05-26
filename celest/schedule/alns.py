

from typing import Any, List, Callable
import copy
import math
import random


_BEST_SCORE = 3
_BETTER_SCORE = 2
_ACCEPT_SCORE = 1
_REJECT_SCORE = 0
_SUM_SCORE = _BEST_SCORE + _BETTER_SCORE + _ACCEPT_SCORE + _REJECT_SCORE


class ALNS:
    """ALNS(initial_solution)

    The ANLS class is a generic framework of the adaptive large
    neighborhood search algorithm.

    Parameters
    ----------
    initial_solution : Any
    """

    def __init__(self, initial_solution: Any) -> None:

        self.initial_solution = initial_solution

    def add_cost_function(self, cost_function: Callable) -> None:
        """Add cost function to minimize.

        Parameters
        ----------
        cost_function : Callable
            Cost function that takes a solution as an input and returns a
            scalar cost value.
        """

        self.cost_function = cost_function

    def add_destroy_functions(self, destroy_functions: List[Callable]) -> None:
        """Add destroy functions to the ALNS instance.

        Parameters
        ----------
        destroy_functions : List[function]
            List of destroy functions that take a solution as an input and
            return a destroyed solution.
        """

        self.destroy_functions = destroy_functions

    def add_repair_functions(self, repair_functions: List[Callable]) -> None:
        """Add repair functions to the ALNS instance.

        Parameters
        ----------
        repair_functions : List[function]
            List of repair functions that take a destroyed solution as an input
            and return a repaired solution.
        """

        self.repair_functions = repair_functions

    def add_completeness_function(self, completeness_function) -> None:
        """Function to check if the solution is complete.

        Parameters
        ----------
        completeness_function : Callable
            Function that takes a solution and returns a boolean.
        """

        self.completeness_function = completeness_function
    
    def solve(self, max_iter: int, p: float, l: float, number_to_remove: int, number_to_add: int) -> Any:
        """Determine the optimal solution.

        Parameters
        ----------
        max_iter : int
            Maximum number of iterations.
        p : float
            Annealing coefficient.
        l : float
            Decay parameter within the range [0, 1].
        number_to_remove : int
            Number of items to remove on each iteration.
        number_to_add : int
            Number of items to add on each iteration.

        Returns
        -------
        Any
            Optimal solution.
        """

        current_solution = copy.deepcopy(self.initial_solution)
        best_solution = copy.deepcopy(self.initial_solution)
        temperature = 0.05 * self.cost_function(current_solution) / math.log(0.5)

        self.destroy_weights = [1] * len(self.destroy_functions)
        self.repair_weights = [1] * len(self.repair_functions)

        if self.completeness_function(best_solution):
            return best_solution

        for _ in range(max_iter):

            i = self._get_destroy_index()
            j = self._get_repair_index()

            temporary_solution = copy.deepcopy(current_solution)
            temporary_solution = self._get_destroy_function(i)(temporary_solution, number_to_remove)
            temporary_solution = self._get_repair_function(j)(temporary_solution, number_to_add)
            score = _REJECT_SCORE

            if self._accept(temporary_solution, current_solution, temperature):
                current_solution = copy.deepcopy(temporary_solution)
                if self.cost_function(temporary_solution) > self.cost_function(current_solution):
                    score = _BETTER_SCORE
                else:
                    score = _ACCEPT_SCORE

            if self.cost_function(temporary_solution) < self.cost_function(best_solution):
                best_solution = copy.deepcopy(temporary_solution)
                score = _BEST_SCORE

            if self.completeness_function(best_solution):
                return best_solution

            self._update_destroy_weights(l, score, i)
            self._update_repair_weights(l, score, j)

            temperature = p * temperature

        return best_solution

    def _get_destroy_index(self) -> int:
        """Return index of the destroy function to be used.

        Returns
        -------
        int
            Index of the destroy function to be used.
        """

        p = self._get_probabilities(self.destroy_weights)
        i = random.choices([j for j in range(len(self.destroy_functions))], p)[0]

        return i
    
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

    def _get_repair_index(self) -> int:
        """Return index of the repair function to be used.

        Returns
        -------
        int
            Index of the repair function to be used.
        """

        p = self._get_probabilities(self.repair_weights)
        i = random.choices([j for j in range(len(self.repair_functions))], p)[0]

        return i

    def _get_destroy_function(self, i) -> Callable:
        """Return destroy function at index i.

        Parameters
        ----------
        i : int
            Destroy function index.

        Returns
        -------
        function
        """

        return self.destroy_functions[i]

    def _get_repair_function(self, i) -> Callable:
        """Return repair function at index i.

        Parameters
        ----------
        i : int
            Repair function index.

        Returns
        -------
        function
        """

        return self.repair_functions[i]
    
    def _accept(self, temporary_solution: Any, current_solution: Any, T: float) -> bool:
        """Accpetance criterion.

        This method implements a simulated annealing acceptance criterion.

        Parameters
        ----------
        temporary_solution : Any
            Temporary solution.
        current_solution : Any
            Current solution.
        T : float
            Temperature.

        Returns
        -------
        bool
            True if the solution is accepted, False otherwise.
        """

        current_cost = self.cost_function(current_solution)
        temporary_cost = self.cost_function(temporary_solution)

        if temporary_cost < current_cost:
            return True
        else:
            probability = math.exp(100 * (temporary_cost - current_cost) / 
                                   (T * current_cost))
            return random.choices([False, True], [1 - probability, probability])[0]

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

        self.destroy_weights[i] = (1 - l) * self.destroy_weights[i] + l * score / _SUM_SCORE

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

        self.repair_weights[i] = (1 - l) * self.repair_weights[i] + l * score / _SUM_SCORE
