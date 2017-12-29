"""Thi module is from Udacity course on Linear Algebra"""
from math import acos, degrees, pi
from decimal import Decimal, getcontext
getcontext().prec = 30


class Vector(object):
    """This class contains the properties
    and modules to create N dimentional vectors"""
    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    NO_UNIQUE_PARALLEL_COMPONENT_MSG = 'No unique parallel component'
    NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG = 'No unique orthogonal component'
    ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG = 'Defined only in two or \
    three dimentions only'
    VECTOR_DIM_MISMATCH_MSG = 'Vectors'' dimentions do not match'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def is_zero(self, tolerance=1e-10):
        """Checks if this vector is zero"""
        return self.magnitude() < tolerance

    def is_parallel_to(self, input_vector):
        """Checks if this vector is parallel to the input vector"""
        return (self.is_zero() or
                input_vector.is_zero() or
                self.angle_with(input_vector) == 0 or
                self.angle_with(input_vector) == pi)

    def is_orthogonal_to(self, input_vector, tolerance=1e-10):
        """Checks if this vector is orthogonal"""
        return abs(self.dot(input_vector)) < tolerance

    def plus(self, input_vector):
        """Sum of this vector with the input"""
        new_coordinates = [x+y for x, y in
                           zip(self.coordinates, input_vector.coordinates)]
        return Vector(new_coordinates)

    def minus(self, input_vector):
        """Subtract the input from the current vector)"""
        new_coordinates = [x-y for x, y in
                           zip(self.coordinates, input_vector.coordinates)]
        return new_coordinates

    def magnitude(self):
        """Get the magnitude of the vector"""
        return sum([x**2 for x in self.coordinates])**Decimal('0.5')

    def normalized(self):
        """Get the normalized form of the vector"""
        try:
            magnitude = self.magnitude()
            return self.times_scalar(Decimal('1.0')/magnitude)
        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def times_scalar(self, multiplier):
        """Get scalar multiplication value of the vector in list format"""
        new_coordinates = [Decimal(multiplier)*x for x in self.coordinates]
        return Vector(new_coordinates)

    def dot(self, input_vector):
        """Get dot product value of this vector with the input vector"""
        return sum((x*y for x, y in
                    zip(self.coordinates, input_vector.coordinates)))

    def angle_with(self, input_vector, in_degrees=False):
        """Get the angle between this vector
        and the input vector in either radian(default) or degrees"""
        vector_one = self.normalized()
        vector_two = input_vector.normalized()
        angle_in_radians = acos(vector_one.dot(vector_two))
        if in_degrees is False:
            return_value = angle_in_radians
        else:
            return_value = degrees(angle_in_radians)
        return return_value

    def component_parallel_to(self, basis):
        """Returns the projection of self vector on the input vector"""
        try:
            normalized = basis.normalized()
            weight = self.dot(normalized)
            return normalized.times_scalar(weight)
        except Exception as exp:
            if str(exp) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise exp

    def component_orthogonal_to(self, basis):
        """Returns the orthogonal vector of itself
        when the input is the basis"""
        try:
            projection = self.component_parallel_to(basis)
            return self.minus(projection)
        except Exception as exp:
            if str(exp) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
            else:
                raise exp

    def cross(self, input_vector):
        """Returns the vector which is the result of the cross product
        of itself and the input vector"""
        vector1 = self.coordinates
        vector2 = input_vector.coordinates
        if len(vector1) > 3 or len(vector2) > 3:
            raise Exception(self.ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG)
        elif len(vector1) != len(vector2):
            raise Exception(self.VECTOR_DIM_MISMATCH_MSG)
        else:
            if len(vector1) == 2:
                vector1 = vector1 + (Decimal('0'),)
                vector2 = vector2 + (Decimal('0'),)
            item0 = (vector1[1]*vector2[2] -
                     vector2[1]*vector1[2])
            item1 = -(vector1[0]*vector2[2] -
                      vector2[0]*vector1[2])
            item2 = (vector1[0]*vector2[1] -
                     vector2[0]*vector1[1])
            return Vector([item0, item1, item2])

    def area_of_parallelogram_with(self, input_vector):
        """Returns the area of the parallelogram created by
        connecting itself with the input vector"""
        cross_result = self.cross(input_vector)
        return cross_result.magnitude()

    def area_of_triangle_with(self, input_vector):
        """Returns the area of the parallelogram created by
        connecting itself with the input vector"""
        return Decimal('0.5')*self.area_of_parallelogram_with(input_vector)
