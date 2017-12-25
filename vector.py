"""Thi module is from Udacity course on Linear Algebra"""
class Vector(object):
    """This class contains the properties and modules to create N dimentional vectors"""
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple(coordinates)
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)
    def __eq__(self, v):
        return self.coordinates == v.coordinates
    def magnitude(self):
        """Get the magnitude of the vector"""
        return sum([x**2 for x in self.coordinates])**0.5
        # coordinates_squared = [x**2 for x in self.coordinates]
        # return sum(coordinates_squared);
    def normalized(self):
        """Get the normalized form of the vector"""
        magnitude = self.magnitude()
        return self.times_scalar(1./magnitude)
    def times_scalar(self, multiplier):
        """Get scalar multiplication value of the vector in list format"""
        return [x*multiplier for x in self.coordinates]
