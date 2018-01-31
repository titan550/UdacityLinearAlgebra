from decimal import Decimal, getcontext
from copy import deepcopy
from vector import Vector
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
                c = MyDecimal(system[row].normal_vector.coordinates[col])
                if c.is_near_zero():
                    swap_with_index = (next((row2 for row2, plane in enumerate(system)
                                            if row2 > row and
                                             not MyDecimal(plane.normal_vector.coordinates[col]).is_near_zero()), None))
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

    def compute_rref(self):
        tf = self.compute_triangular_form()
        num_equations = len(tf)
        pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()
        for i in range(num_equations)[::-1]:
            j = pivot_indices[i]
            if j > 0:
                tf.scale_row_to_make_coefficient_equal_one(i, j)
                tf.clear_coefficients_above(i, j)
        return tf

    def scale_row_to_make_coefficient_equal_one(self, row, col):
        coefficient = Decimal('1.0') / self[row].normal_vector.coordinates[col]
        self.multiply_coefficient_and_row(coefficient, row)

    def clear_coefficients_above(self, row, col):
        for k in range(row)[::-1]:
            alpha = - self[k].normal_vector.coordinates[col]
            # The above assumes that the coefficient of the pivot variable is = 1, otherwise, we would divide alpha by
            # the coefficient
            self.add_multiple_times_row_to_row(alpha, row, k)

    def compute_solution(self):
        try:
            return self.do_gaussian_elimination_and_extract_solution()
        except Exception as e:
            if(str(e) == self.NO_SOLUTIONS_MSG or
                    str(e) == self.INF_SOLUTIONS_MSG):
                return str(e)
            else:
                raise e

    def do_gaussian_elimination_and_extract_solution(self):
        rref = self.compute_rref()

        rref.raise_exception_if_contradictory_equation()
        rref.raise_exception_if_too_few_pivots()

        num_variables = rref.dimension
        solution_coordinates = [rref.planes[i].constant_term for i in range(num_variables)]
        return Vector(solution_coordinates)

    def raise_exception_if_contradictory_equation(self):
        for p in self.planes:
            try:
                p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == 'No nonzero elements found':
                    constant_term = MyDecimal(p.constant_term)
                    if not constant_term.is_near_zero():
                        raise Exception(self.NO_SOLUTIONS_MSG)
                else:
                    raise e

    def raise_exception_if_too_few_pivots(self):
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index >= 0 else 0 for index in pivot_indices])
        num_variables = self.dimension

        if num_pivots < num_variables:
            raise Exception(self.INF_SOLUTIONS_MSG)

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
