"""This module is from the Udacity Course on Linear Algebra"""
from decimal import Decimal, getcontext

from vector import Vector

getcontext().prec = 30


class Line(object):
    """This class defines a line"""
    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = len(normal_vector.coordinates)

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()

    def set_basepoint(self):
        """Finds the basepoint for the line"""
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension

            initial_index = Line.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e

    def intersection_with(self, line):
        if self.dimension != line.dimension or self.dimension != 2:
            raise Exception(NotImplemented)
        if self.is_parallel_to(line):
            if self == line:
                return self
            else:
                return None
        denom = (self.normal_vector[0]*line.normal_vector[1] -
                 self.normal_vector[1]*line.normal_vector[0])
        x_1 = (line.normal_vector[1]*self.constant_term -
               self.normal_vector[1]*line.constant_term)/denom
        x_2 = (self.normal_vector[0]*line.constant_term -
               line.normal_vector[0]*self.constant_term)/denom
        return [x_1, x_2]

    def is_parallel_to(self, line):
        """Returns True if this line is parallel with the input line"""
        return self.normal_vector.is_parallel_to(line.normal_vector)

    def __eq__(self, line):
        """Returns True is the this line is equal to the input line"""
        if self.normal_vector.is_zero():
            if not line.normal_vector.is_zero():
                return False
            else:
                diff = self.constant_term - line.constant_term
                return MyDecimal(diff).is_near_zero()
        elif line.normal_vector.is_zero():
            return False
        if not self.is_parallel_to(line):
            return False
        basepints_vector = self.basepoint.minus(line.basepoint)
        return basepints_vector.is_orthogonal_to(self.normal_vector)

    def __str__(self):
        """Creates the string format of this line's equation"""
        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = Line.first_nonzero_index(n)
            terms = [write_coefficient(n[i],
                                       is_initial_term=(i == initial_index)) +
                     'x_{}'.format(i+1)
                     for i in range(self.dimension)
                     if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output

    @staticmethod
    def first_nonzero_index(iterable):
        """Returns the index of the first non-zero
        element in the input vector"""
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)


class MyDecimal(Decimal):
    """Extension class for decimal operations"""
    def is_near_zero(self, eps=1e-10):
        """Checks if itself is smaller than the
        input value or the default value"""
        return abs(self) < eps
