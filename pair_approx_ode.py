'''
This files contains functions needed to compute the mean-field and pair-approximation for the paper

The power of two choices on graphs: the pair-approximation is accurate?
N Gast - ACM SIGMETRICS Performance Evaluation Review, 2015 

Main functions:
- mm1(rho, N)

'''
import numpy as np

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

print(mm1(0.7))

def twoChoiceTheory(rho, d=50):
    """
    Returns the steady-state distribution of the mean field approx. of
    a "power of two choice" model.
    """
    assert rho>=0 and rho<1, "rho should be in [0,1)"
    x = rho ** (2**np.arange(d)-1)
    return -np.diff(x)

print(twoChoiceTheory(0.9))

def derivative_2neighboorChoice(y, rho, mu, N=50):
    """
    Returns the derivative for the pair-approximation.

    Here y[i,j]=percentage of queue i with neighboor j (as defined as in the paper).
    """
    derivative = np.zeros( (N, N))
    p = np.zeros(N)

    for i in range(N):
        if sum(y[i]) >0:
            p[i] = (sum (y[i][i+1:N]) + y[i][i]/2) / sum(y[i])
    for i in range(N):
        for j in range(N):
            if i>0:
                derivative[i][j]   +=  rho*p[i-1]*y[i-1][j] - mu*y[i][j]
                derivative[i-1][j] += -rho*p[i-1]*y[i-1][j] + mu*y[i][j]
                if (i<=j):
                    derivative[i][j]   +=  rho*y[i-1][j]
                    derivative[i-1][j] += -rho*y[i-1][j]
                elif(i-1==j):
                    derivative[i][j]   +=  rho*y[i-1][j]/2
                    derivative[i-1][j] += -rho*y[i-1][j]/2
            if j>0:
                derivative[i][j]   +=  rho*p[j-1]*y[i][j-1] - mu*y[i][j]
                derivative[i][j-1] += -rho*p[j-1]*y[i][j-1] + mu*y[i][j]
                if (j<=i):
                    derivative[i][j]   +=  rho*y[i][j-1]
                    derivative[i][j-1] += -rho*y[i][j-1]
                elif (i==j-1):
                    derivative[i][j]   +=  rho*y[i][j-1]/2
                    derivative[i][j-1] += -rho*y[i][j-1]/2
    return derivative



def stationnary_distribution_pairApprox(rho, N):
    filename = ''.join(['results/pair_approx_N',str(N),'_',str(rho)])
    try:
        x = loadtxt(filename)
        print('we load the saved result')
        return x[:,1]
    except Exception as e:
        print(e)
        print('computation was not done before. We do it and save the results')
    
    y = array([[0.]*N]*N)
    for i in range(0,10):
        for j in range(0,10):
            y[i][j]  = (rho**i)*(rho**j)
    y = y/sum(y)
    h = 0.1
    derivative = array([[0.]*N]*N)
    x0 = -1
    x5 = -1
    T = 100
    conv = 1
    while (abs(sum(y[0])-1+rho)>1e-5
           or abs(x0-sum(y[0])) > 1e-5
           or abs(x5-sum(y[5])) > 1e-7):
        #print( abs(sum(y[0])-1+rho))
        x0 = sum(y[0])
        x5 = sum(y[5])
        for t in range(0,T):
            derivative_2neighboorChoice(y,rho,1,N, derivative)
            y  = y  + h * derivative
        print('T=',T*conv,'x0=',x0,'x5=',x5)
        conv += 1
        # for i in range(0,10):
        #     x[t,i] = sum(y[i])
        #x[t] = sum(y[0])
    x = array([[0.,0.]]*N)
    for i in range(0,N):
        x[i,0] = i
        x[i,1] = sum(y[i])
    savetxt(filename, x)
    return (x[:,1])

