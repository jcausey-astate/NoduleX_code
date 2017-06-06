import numpy as np

def bounding_cube(arr):
    '''
    Report the limits of the smallest box that contains the all values that
    boolean-evaluate to True. Returns a series of three tuples, corresponding
    to the limits along each of the axes 0, 1, and 2.
    '''
    assert len(arr.shape) == 3
    a0_columns = np.any(arr, axis=0)
    (a2lim, a1lim) = bounding_box(a0_columns)

    (a2lim_extra, a0lim) = bounding_box(np.any(arr, axis=1))

    assert a2lim == a2lim_extra
    assert bounding_box(np.any(arr, axis=2)) == (a1lim, a0lim)

    return (a0lim, a1lim, a2lim)


def bounding_box(arr):
    '''
    Determine the smallest box that contains all the true-valued elements
    in the given two dimensional matrix.
    '''
    assert len(arr.shape) == 2

    x_line = np.any(arr, axis=0)
    xlim = bounding_line(x_line)

    y_line = np.any(arr, axis=1)
    ylim = bounding_line(y_line)

    return (xlim, ylim)


def bounding_line(arr):
    '''
    Determine the first and last true value in the given boolean-valued
    vector
    '''
    assert len(arr.shape) == 1

    first_true = np.argmax(arr)
    last_true = len(arr) - np.argmax(arr[::-1])

    return (first_true, last_true)
