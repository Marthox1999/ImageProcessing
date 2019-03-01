import numpy as np

def get_gaussian_filter(neighbours = 1, sigma = 1):
    '''
    This function calculates a NxN gaussian filter; N = neighbours * 2 + 1
    ----------
    G(X,Y) = (1/ 2 * PI * SIGMA) * e^(-(X^2 + Y^2) / 2 * SIGMA^2)
    ----------
    The values of X and Y are the difference between the current pixel and the center of the kernel
    For a 3x3 kernel, we have the following values of X and Y paired as tuples.
    ----------
    [
        [(-1,-1), (0, -1) (1, -1)],
        [(-1, 0), (0, 0), (1, 0)],
        [(-1, 1), (0, 1), (1, 1)]
    ]
    '''
    N =  neighbours * 2 + 1

    #First, we create an array of shape NxN
    X = Y = np.zeros((N, N))

    #Then, we populate all the columns using np.arange - neighbours, to get the distance
    #between the element in the middle (index = neigboors) and the index given by arange

    #NOTE: If you want to change the indexes, similar to the Rayleigh function, change to the following line:
    #Y[:] = np.arange(N, dtype=np.float32)
    Y[:] = np.arange(N, dtype=np.float32) - neighbours
    
    #To calculate the rows, we transpose the columns
    X = np.transpose(Y)

    #NOTE: We do this to use numpy's vectorization, since we can have X and Y index values separated. 
    #Remember that the Gaussian filter in 2D is equal to multiplying two 1D Gaussians

    gaussian_filter = np.zeros((N, N))

    #We write down the left side (scalar value) of the equation: (1/ 2 * PI * SIGMA)
    left_side = (1 / (np.pi * 2 * np.power(sigma, 2))) 

    #We start writting down the right side: -(X^2 + Y^2) / 2 * SIGMA^2
    aux = -(np.power(X, 2) + np.power(Y, 2)) / (2 * np.power(sigma, 2))
    #Finally, we use the Exp function to calculate the e^(aux)
    right_side = np.exp(aux)

    #Multiply both sides of the equation
    gaussian_filter = left_side * right_side

    #To calculate the scalar value to multiply by, we sum over all the axes of the gaussian matrix, using np.sum
    factor_value = 1 / np.sum(gaussian_filter)

    return gaussian_filter, factor_value

def get_rayleigh_filter(neighbours = 1, sigma = 1):
    '''
    This function calculates a NxN rayleigh filter; N = neighbours * 2 + 1
    ----------
    G(X,Y) = [X * Y / SIGMA^4] * [e^(-(X^2 + Y^2) / 2 * SIGMA^2))]
    ----------
    The values of X and Y are the indixes of each pixel
    For a 3x3 kernel, we have the following values of X and Y paired as tuples.
    ----------
    [
        [(0,0), (0, 1) (0, 2)],
        [(1, 0), (1, 1), (1, 2)],
        [(2, 0), (2, 1), (2, 2)]
    ]
    '''
    N =  neighbours * 2 + 1

    #First, we create an array of shape NxN
    X = Y = np.zeros((N, N))

    #Then, we populate all the columns using np.arange to get the indexes
    Y[:] = np.arange(N, dtype=np.float32)
    
    #To calculate the rows, we transpose the columns
    X = np.transpose(Y)

    #NOTE: We do this to use numpy's vectorization, since we can have X and Y index values separated. 
    #Similar to gaussian 2D Rayleighs are equal to multiplying two 1D Rayleighs
    rayleigh_filter = np.zeros((N, N))

    #We write down the left side of the equation: [X * Y / SIGMA^4]
    left_side = np.divide(X * Y, np.power(sigma, 4))

    #We start writting down the right side: -(X^2 + Y^2) / 2 * SIGMA^2
    aux = -(np.power(X, 2) + np.power(Y, 2)) / (2 * np.power(sigma, 2))
    #Finally, we use the Exp function to calculate the e^(aux)
    right_side = np.exp(aux)


    #Multiply both sides of the equation
    rayleigh_filter = left_side * right_side

    #To calculate the scalar value to multiply by, we sum over all the axes of the gaussian matrix, using np.sum
    factor_value = 1 / np.sum(rayleigh_filter)

    return rayleigh_filter, factor_value

def get_integer_valued_gaussian_filter(neighbours = 1, sigma = 1):
    '''
    This function approximates to a gaussian distribution using the N-th pascal triangle
    '''
    N = neighbours * 2 + 1

    #First, we calculate the N-th row of the pascal triangle. The function is exlcusive, so we substract by 1
    #The results is the converted to a numpy's ndarray
    pascal_row = np.asarray(get_kth_pascal_row(N - 1), dtype=np.int32)
    #We expand the dimension of the row, this returns the transpose of the row
    pascal_row = np.expand_dims(pascal_row, 1)
    
    #We create the extended version by replicating the rows len(row) times
    #NOTE: See https://dsp.stackexchange.com/questions/10057/gaussian-blur-standard-deviation-radius-and-kernel-size
    extended_pascal = np.resize(pascal_row, (pascal_row.shape[0], pascal_row.shape[0]))    

    #We perform a pointwise multp to calculate our result
    result = pascal_row * extended_pascal
    #We calculate the scalar factor, which is the sum of all of the values in the matrix
    scalar_factor = np.sum(result)

    return result, scalar_factor

def get_kth_pascal_row(row_number):
    '''
    Parameters

    @row_number: the row number to return, it is exclusive (starts at 0)
    ---------
    This function calculates the k-th row of the pascal triangle.
    '''
    if row_number == 0:
        return [1, ]
        
    last_row = [1, ]
    for R in range(1, row_number+1):
        row = []
        row.append(1)
        for C in range(R - 1):
            row.append(last_row[C] + last_row[C+1])
        row.append(1)
        
        last_row = row
    return last_row
