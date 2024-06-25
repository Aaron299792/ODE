def Euler(x0,t,function):
    
    """
    Solves ordinary differential equations (ODEs) using Euler's method for a given initial condition for a defined number of time steps (independent variable steps),for each timestep pass as argument in an array, the function computes the time step size, then it allocates and iterates over an array using Euler's method, such array is then return by the function. Its important to take into account that for the Euler's method the computation error is liniar in the size of the time step, so the closed the time steps are the better results are obtained.

    -it pressumes the use of a numpy array in the arguments so it is necessary to import the numpy module for its use.
    -if t is not an rank 1 numpy array, the function will not work.
    -time steps refers to any independent variable steps in these case.

    Examples:
        >>> Euler(1.0, numpy.array[0.0, 0.1, 0.2, 0.3, 0.4], cuadratic_function)
        array([1.        , 1.1       , 1.221     , 1.3700841 , 1.55779714])

    Args:
    x0 (float): First argument, initial condition of a given ODE problem.
    t  (numpy array): Second argument, 1D numpy array that contains the time steps to be evaluated in the solution.
    function : (block statement that returns a float): Third argument, function that evaluates the 1st order derivative for two given values.

    Returns:
    output (numpy array):
        Returns an 1D numpy array containing the dependent variable evaluated in a given time step using the Euler's method

    """

    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0
    for i in range(t.size - 1):
        x[i+1] = x[i] + (h * f(x[i], t[i]))
    return x

def Rk2(x0, t, f):

    """ 
    Solves ordinary differential equations (ODEs) using a Runge-Kutta method of 2nd order (RK2)for a given initial condition for a defined number of time steps (independent variable steps),for each timestep pass as argument in an array, the function computes the time step size, then it allocates and iterates over an array using the method, such array is then return by the function. Its important to take into account that for the Runge-Kutta of 2nd order the more time steps the more precise results are get.

    -it pressumes the use of a numpy array in the arguments so it is necessary to import the numpy module for its use.
    -if t is not an rank 1 numpy array, the function will not work.it pressumes the use of a numpy array in the arguments
    -time steps refers to any independent variable steps in these case.

    Examples:
        >>> Rk2( 1.0, numpy.array[0.0, 0.1, 0.2, 0.3, 0.4] , cuadratic_function)
        array([1.        , 0.91926042, 0.86479208, 0.82975512, 0.80971762])

    Args: 
    x0 (float): First argument, initial condition of a given ODE problem.
    t (array_like): Second argument, 1D numpy array that contains the time steps to be evaluated in the solution of the ODE.
    function (block statement that returns a float): Third argument, function that evaluates the 1st order derivative for two given values.

    Returns
    output (array_like):  
        Returns an 1D numpy array containing the dependent variable evaluated in a given time step using the Runge-Kutta method of 2nd order.

    """

    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0

    for i in range (t.size-1):
        k1 = h*func(x[i],t[i])
        k2 = h*func(x[i] + 0.5*k1 , t[i] + 0.5*h)
        x[i+1] = x[i] + k2

    return x

def Rk4(x0,t,f):

    """ 
    Solves ordinary differential equations (ODEs) using a Runge-Kutta method 4th (RK4)for a given initial condition for a defined number of time steps (independent variable steps),for each timestep pass as argument in an array, the function computes the time step size, then it allocates and iterates over an array using the Runge-Kutta method of 4th, such array is then return by the function. Its important to take into account that for this method the more time steps and shorter time step size, the better the results are

    -it pressumes the use of a numpy array in the arguments so it is necessary to import the numpy module for its use.
    -if t is not an rank 1 numpy array, the function will not work.it pressumes the use of a numpy array in the arguments
    -time steps refers to any independent variable steps in these case.
    
    Examples:
        >>> Rk4(1.0, numpy.array[0.0, 0.1, 0.2, 0.3, 0.4], cuadratic_function)
        array([1.        , 1.11111049, 1.24999799, 1.42856619, 1.66665326])

    Args:
    x0 (float): First argument, initial condition of a given ODE problem.
    t (array_like): Second argument, 1D numpy array that contains the time steps to be evaluated in the solution of the ODE.
    function (block statement that returns a float): Third argument, function that evaluates the 1st order derivative for two given values.

    Returns
    output (array_like):
        Returns an 1D numpy array containing the dependent variable evaluated in a given time step using the Runge-Kutta method of 4th order.
    
    """
    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0

    for i in range (t.size-1):
        k1 = h*f(x[i],t[i])
        k2 = h*f(x[i] + 0.5*k1, t[i] + 0.5*h)
        k3 = h*f(x[i] + 0.5*k2, t[i] + 0.5*h)
        k4 = h*f(x[i] + k3, t[i] + h)
        x[i+1] =  x[i] + (1/6)*(k1 + 2.0*k2 + 2.0*k3 + k4)
    return x
