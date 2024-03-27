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

def Gaussian_filter(x, y, sigma, normalize = True, normalize_to = 1.0):
    """Use convolution to smooth the data y(x) with a standard deviation sigma."""
    # check if uniform grid of x, otherwise raise an error
    dx = x[1] - x[0]
    if not np.allclose(np.diff(x), dx):
        raise ValueError("The x data is not uniformly spaced.")
    assert len(x) == len(y)
    # if full of y < 0, first reflect the data
    reflected = False
    if np.all(y <= 0):
        y = -y
        reflected = True

    y_smoothed = sp.ndimage.filters.gaussian_filter1d(y, sigma/dx)
    if normalize:
        norm = sp.integrate.simps(y_smoothed, x)
        y_smoothed = y_smoothed * normalize_to / norm

    y_smoothed = y_smoothed if not reflected else -y_smoothed

    return y_smoothed