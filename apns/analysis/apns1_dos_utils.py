import numpy as np
import scipy as sp

def zero_padding(xmin: float, xmax: float, dx: float, x: np.ndarray, y: np.ndarray):
    """Zero padding for the x and y data, so that the x data is within the range [xmin, xmax] with a step of dx.
    The y data is interpolated to the new x data.
    
    Args:
        xmin (float): minimum value of the new x data
        xmax (float): maximum value of the new x data, EXCLUSIVE!
        dx (float): step of the new x data
        x (np.ndarray): original x data
        y (np.ndarray): original y data
        interpolate (str): interpolation method, "none", "linear", "cubic"
    
    Returns:
        np.ndarray: new x data
        np.ndarray: new y data
    """
    xnew = np.arange(xmin, xmax, dx)
    ynew = np.zeros(len(xnew))

    def nearest_index(x, xnew):
        """return the nearest index of x in xnew, and the difference between x and xnew[index]"""
        i = np.abs(xnew - x).argmin()
        return i, xnew[i] - x
    
    for i in range(len(xnew)):
        if xnew[i] < x[0] or xnew[i] > x[-1]:
            continue
        _i, d = nearest_index(xnew[i], x)
        if abs(d) < dx:
            ynew[i] = y[_i]
    return xnew, ynew

import unittest
class APNS1DOSUtilitiesTest(unittest.TestCase):
    def test_zero_padding(self):
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 2, 3, 4, 5])
        xnew, ynew = zero_padding(0, 7, 1, x, y)
        self.assertTrue(np.allclose(xnew, np.array([0, 1, 2, 3, 4, 5, 6])))
        self.assertTrue(np.allclose(ynew, np.array([0, 1, 2, 3, 4, 5, 0])))

        x = np.linspace(0, 2 * np.pi, 5)
        y = np.sin(x)
        xnew, ynew = zero_padding(0, 2 * np.pi, np.pi / 10, x, y)
        self.assertTrue(np.allclose(xnew, np.linspace(0, 2 * np.pi, 20, endpoint=False)))
        yref = [0.0000000e+00, 0.0000000e+00, 0.0000000e+00,  0.0000000e+00,
                0.0000000e+00, 1.0000000e+00, 0.0000000e+00,  0.0000000e+00,
                0.0000000e+00, 0.0000000e+00, 1.2246468e-16,  0.0000000e+00,
                0.0000000e+00, 0.0000000e+00, 0.0000000e+00, -1.0000000e+00,
                0.0000000e+00, 0.0000000e+00, 0.0000000e+00,  0.0000000e+00]
        self.assertEqual(len(ynew), len(yref))
        self.assertTrue(np.allclose(ynew, yref))

if __name__ == "__main__":
    unittest.main()