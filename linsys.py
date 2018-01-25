from decimal import Decimal, getcontext
from copy import deepcopy
# from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live\
    in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def swap_rows(self, row1, row2):
        self.planes[row1], self.planes[row2] = (self.planes[row2],
                                                self.planes[row1])

    def multiply_coefficient_and_row(self, coefficient, row):
        normal_vector = self[row].normal_vector.times_scalar(coefficient)
        constant_term = self[row].constant_term * coefficient
        self[row] = Plane(normal_vector, constant_term)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add,
                                      row_to_be_added_to):
        n1 = self[row_to_add].normal_vector
        n2 = self[row_to_be_added_to].normal_vector
        k1 = self[row_to_add].constant_term
        k2 = self[row_to_be_added_to].constant_term
        new_normal_vector = n1.times_scalar(coefficient).plus(n2)
        new_constant_term = (k1 * coefficient) + k2
        self[row_to_be_added_to] = Plane(new_normal_vector, new_constant_term)

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations
        for i, p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def compute_triangular_form(self):
        system = deepcopy(self)
        col = 0
        for row in range(len(system)):
            while col < system.dimension:
                if system[row].normal_vector.coordinates[col] == 0:
                    swap_with_index = (next((row2 for row2, plane in enumerate(system)
                                            if row2 > row and plane.normal_vector.coordinates[col] != 0), None))
                    # Gets the first plane with a non-zero value at dimension i
                    if swap_with_index:
                        system.swap_rows(row, swap_with_index)
                    else:
                        col += 1
                        continue
                for row2 in range(row+1, len(system)):
                    coefficient = - system[row2].normal_vector[col]/system[row].normal_vector[col]
                    system.add_multiple_times_row_to_row(coefficient, row, row2)
                # Sets all the following rows to zero
                col += 1
                break
        return system

    def __len__(self):
        return len(self.planes)

    def __getitem__(self, i):
        return self.planes[i]

    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1, p) for i, p in
                enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps
