# Python ODE solvers implementation

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # Introduction to the numerical methods to solve ODEs. 
        reference.md  # The documentation of the methods implemented in python.

## Numerical methods to solve ordinary differential equations (ODEs)

### Euler's method

This method is derived directly from the Taylor's first order expansion of a function \(x(t)\). Assuming a finite increment \(h\), we get

$$
x(t + h) = x(t) + h\frac{dx}{dt} + \frac{h^2}{2}\frac{d^2x}{dt^2} + O(h^3)
$$

The case being \(h\) sufficently small, we can get the function by means of 

$$
x(t + h) = x(t) + hf(x,t)
$$

where $$f(x,t) = \frac{dx}{dt}$$

The Euler's algorithim works as follows:

* We star off with \(t = t_0\),  \(x = x_0\), the first one being the initial time and the second one the function we seek evaluated at that time.
* An structure that contains a discretization of the time is created for equally distance time steps, such distance is \(h\).
* For each time in our structure we find a x using the discretize form of the previous equation:
$$x_i = x_{i-1} + hf(x_{i-1}, t_{i-1})$$

For this method is easy from the Taylor expansion to se that the error grows as $$\frac{h}{2}[f_b-f_a]$$

### Runge-Kutta methods
This is a family of methods of different order that uses the Taylor's first order expansion taking into acount more than one previous step to compute the first order differential equation this allow to improve the approximation. In comparison to using only one step behind this family of methods converges faster to the solution than the Euler's method. 

#### Runge-Kutta method of 2nd order (RK2)

This method uses an approximation in half size steps compare to the full step that uses the Euler method such that we evaluate the function in \(t = t + h/2\). This improves the aproximation for the same value of h that is used in the Euler method. 

Approximating the function through the first order term of the Taylor's expansion we get:

$$
x(t + h) = x(t) + hf[x(t+h/2), t + h/2] + O(h^3)
$$

Finally we implement the following equations in our algorithim to implement the RK2 method:


* $$k_1 = hf(x,t)$$
* $$k_2 = hf\left(x + \frac{k_1}{2}, t + \frac{h}{2}\right)$$ 
* $$ x(t + h) = x(t) + k_2 $$

#### Runge-Kutta method of 4th order (RK4)
Using the same approach as the Euler's and RK2 method, we find a expresion for the  function using an expansion of the Taylor's series using 4 steps rather than one, this allows after a cumbersome algebra to get ride of the \(h^3\) and \(h^4\) orders and to get a approximation with an error that grows as \(h^5\), this makes the RK4 a method which is consider for many the perfect balance between complexity and error of approximation. The RK4 equations are:

* $$k_1 = hf(x,t)$$
* $$k_2 = hf(x+\frac{k_1}{2}, t + \frac{h}{2})$$
* $k_3 = hf(x + \frac{k_2}{2}, t + \frac{h}{2})$$
* $$k_4 = hf(x + k_3, t +h)$$
* $$x(t + h) = x(t) + \frac{1}{6}(k_1  + 2k_2 + 2k_3 + k_4)$$

This method is useful in most applications, it's easy to code and the results are precisse. 
