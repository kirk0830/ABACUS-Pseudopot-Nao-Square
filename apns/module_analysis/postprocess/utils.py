import numpy as np

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

def Gauss_smoothing(x, y, sigma, normalize = True, normalize_to = 1.0):
    """Gaussian smoothing of the data y(x) with a standard deviation sigma.

    Args:
        x (np.ndarray): x-axis values.
        y (np.ndarray): y-axis values.
        sigma (float): standard deviation of the Gaussian.

    Returns:
        np.ndarray: smoothed y-axis values.
    """
    n = len(x)
    m = len(y)
    assert n == m
    y_smoothed = np.zeros(n)
    for i in range(n):
        for j in range(n):
            y_smoothed[i] += y[j]*np.exp(-(x[i] - x[j])**2/(2*sigma**2))
    if normalize:
        norm = np.trapz(y_smoothed, x)
        y_smoothed = y_smoothed * normalize_to / norm
            
    return y_smoothed
