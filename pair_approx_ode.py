'''
This files contains functions needed to compute the mean-field and pair-approximation for the paper

The power of two choices on graphs: the pair-approximation is accurate?
N Gast - ACM SIGMETRICS Performance Evaluation Review, 2015 

Main functions:
- mm1(rho, d=50)
- twoChoiceTheory(rho, d=50)
- 
'''
import numpy as np
import scipy.integrate
try:
    from numba import jit
except ImportError:
    print("Numba is not installed. We will run without it")
    print("This can be *much* slower")
    def jit(nopython=True):
        def useless_decorator(func):
            return func
        return useless_decorator

class NEIGHBOOR_TWO_CHOICE():
    """
    proxy class to call various functions to compare simulation and approximations.
    """
    def __init__(self, rho, d=50):
        """
        d (default=50) should not a priori be changed.
        """
        self.rho = rho
        self.d = d

    def mm1(self):
        """
        Return the no choice value
        """
        return mm1(self.rho, self.d)

    def mean_field_approximation(self):
        """
        Returns the mean field approximation
        """
        return twoChoiceTheory(self.rho, self.d)

    def simulation(self):
        assert False, "not implemented"
        pass

    def pair_approximation(self, force_recompute=False, atol=1e-10):
        """
        Returns the pair approximation. This methods numerically integrate
        an ODE corresponding to the fixed point.

        Args:
        - force_recompute: if False, will try to load a file from disk
        - atol: absolute numerical error tolerated (when computing the fixed point of the ODE)
        """
        return stationnary_distribution_pairApprox(self.rho, self.d, force_recompute, atol)


def mm1(rho, d=50):
    """
    Returns the steady-state distribution of a M/M/1.
    
    Output:
    - an array of size N where x[i] = proba(queue=i job)

    Args:
    - rho = load of the M/M/1.
    - d (=50): limit of the queue lengths
    """
    assert rho>=0 and rho<1, "rho should be in [0,1)"
    x = rho**np.arange(d+1)
    return -np.diff(x)

def twoChoiceTheory(rho, d=50):
    """
    Returns the steady-state distribution of the mean field approx. of
    a "power of two choice" model.
    """
    assert rho>=0 and rho<1, "rho should be in [0,1)"
    x = rho ** (2**np.arange(d)-1)
    return -np.diff(x)

@jit(nopython=True)
def derivative_2neighboorChoice(y, rho, mu=1, d=50):
    """
    Returns the derivative for the pair-approximation.

    Here y[i,j]=percentage of queue i with neighboor j (as defined as in the paper).
    """
    derivative = np.zeros( (d, d))
    p = np.zeros(d)

    for i in range(d):
        if np.sum(y[i]) >0:
            p[i] = (np.sum (y[i][i+1:d]) + y[i][i]/2) / np.sum(y[i])
    for i in range(d):
        for j in range(d):
            if i>0:
                derivative[i][j]   +=  rho*p[i-1]*y[i-1][j] - mu*y[i][j]
                derivative[i-1][j] += -rho*p[i-1]*y[i-1][j] + mu*y[i][j]
                if  i<=j:
                    derivative[i][j]   +=  rho*y[i-1][j]
                    derivative[i-1][j] += -rho*y[i-1][j]
                elif i-1==j:
                    derivative[i][j]   +=  rho*y[i-1][j]/2
                    derivative[i-1][j] += -rho*y[i-1][j]/2
            if j>0:
                derivative[i][j]   +=  rho*p[j-1]*y[i][j-1] - mu*y[i][j]
                derivative[i][j-1] += -rho*p[j-1]*y[i][j-1] + mu*y[i][j]
                if j<=i:
                    derivative[i][j]   +=  rho*y[i][j-1]
                    derivative[i][j-1] += -rho*y[i][j-1]
                elif i==j-1:
                    derivative[i][j]   +=  rho*y[i][j-1]/2
                    derivative[i][j-1] += -rho*y[i][j-1]/2
    return derivative


def stationnary_distribution_pairApprox(rho, d=50, force_recompute=False, atol=1e-10):
    """
    Returns the stationary distribution for the pair-approximation.

    This functin is based on a numerical integration of the ODE.
    """
    filename = ''.join(['results/pair_approx_d',str(d),'_',str(rho)])
    if not force_recompute:
        try:
            return np.loadtxt(filename)
        except Exception as excep:
            print(excep)
            print('computation was not done before. We do it and save the results')

    # We first initialize the ODE
    y = np.zeros((d,d))
    for i in range(0,d):
        for j in range(0,d):
            y[i][j]  = (rho**i)*(rho**j)
    y = y/np.sum(y)
    h = 0.1 # Step of the ODE
    old_y = np.zeros((d,d))
    T = 100
    conv = 1
    while abs(sum(y[0])-1+rho)>1e-5 or np.sum(np.abs(y-old_y)) > atol:
        old_y = np.copy(y)
        for t in range(0,T):
            y  += h * derivative_2neighboorChoice(y, rho, 1, d)
        print('T=', T*conv, 'diff=', np.sum(np.abs(y-old_y)))
        conv += 1
    x = np.sum(y, 1)
    np.savetxt(filename, x)
    return x
