from pylab import *;
#   no choice:

# transitions are:
#  For departure:
#   (i,j) becomes (i-1,j) at rate mu.
#   (i,j) becomes (i,j-1) at rate mu.
#  For arrival:
#   (i,j) becomes (i,j+1) at rate lambda
#   (i,j) becomes (i+1,j) at rate lambda

# ODE no choice:

N = 50; 

def mm1(lam,N):
    y = lam**array(range(0,N));
    x = array([[0.,0.]]*N);
    for i in range(0,N):
        x[i,0] = i;
        x[i,1] = y[i]/sum(y);
    return(x);

def derivative_noChoice(y,lam,mu, N, derivative):
    for i in range(0,N):
        for j in range(0,N):
            derivative[i][j] = 0;
    for i in range(0,N):
        for j in range(0,N):
            if (i>0):
                derivative[i][j]   +=  1*(lam*y[i-1][j] - mu*y[i][j]);
                derivative[i-1][j] += 1*(-lam*y[i-1][j] + mu*y[i][j]);
            if (j>0):
                derivative[i][j]   +=   1*(lam*y[i][j-1] - mu*y[i][j]);
                derivative[i][j-1] += 1*(-lam*y[i][j-1] + mu*y[i][j]);

def derivative_2neighboorChoice(y,lam,mu, N, derivative):
    for i in range(0,N):
        for j in range(0,N):
            derivative[i][j] = 0;
    p = array([0.0]*N);
    
    for i in range(0,N):
        xi = sum(y[i]);
        if (xi>0):
            p[i] = (sum (y[i][i+1:N]) + y[i][i]/2) / xi;    
    for i in range(0,N):
        for j in range(0,N):
            if (i>0):
                derivative[i][j]   +=  lam*p[i-1]*y[i-1][j] - mu*y[i][j];
                derivative[i-1][j] += -lam*p[i-1]*y[i-1][j] + mu*y[i][j];
                if (i<=j):
                    derivative[i][j]   +=  lam*y[i-1][j];
                    derivative[i-1][j] += -lam*y[i-1][j];
                elif(i-1==j):
                    derivative[i][j]   +=  lam*y[i-1][j]/2;
                    derivative[i-1][j] += -lam*y[i-1][j]/2;
            if (j>0):
                derivative[i][j]   +=  lam*p[j-1]*y[i][j-1] - mu*y[i][j];
                derivative[i][j-1] += -lam*p[j-1]*y[i][j-1] + mu*y[i][j];
                if (j<=i):
                    derivative[i][j]   +=  lam*y[i][j-1];
                    derivative[i][j-1] += -lam*y[i][j-1];
                elif (i==j-1):
                    derivative[i][j]   +=  lam*y[i][j-1]/2;
                    derivative[i][j-1] += -lam*y[i][j-1]/2;

def twoChoiceTheory(lam,N):
    # we assume mu=1
    x = [0]*(N+1);
    rho = lam;
    for i in range(0,N+1):
        x[i] = rho/lam;
        rho = rho*rho;
    y = array([[0.,0.]]*N);
    for i in range(0,N):
        y[i,0] = i;
        y[i,1] = x[i]-x[i+1];
    return y; #(-diff(x));

#y = array([[0.]*N]*N);  y[0][0]  = 1.;
#y2 = array([[0.]*N]*N); y2[0][0] = 1.;

#lam = .7;
#mu = 1;


def stationnary_distribution_pairApprox(lam,N):
    filename = ''.join(['results/pair_approx_N',str(N),'_',str(lam)]);
    try:
        x = loadtxt(filename);
        print('we load the saved result');
        return x[:,1];
    except Exception as e:
        print(e);
        print('computation was not done before. We do it and save the results');
    
    y = array([[0.]*N]*N);
    for i in range(0,10):
        for j in range(0,10):
            y[i][j]  = (lam**i)*(lam**j);
    y = y/sum(y);
    h = 0.1;
    derivative = array([[0.]*N]*N);
    x0 = -1; x5 = -1;
    T = 100;
    conv = 1;
    while (abs(sum(y[0])-1+lam)>1e-5
           or abs(x0-sum(y[0])) > 1e-5
           or abs(x5-sum(y[5])) > 1e-7):
        #print( abs(sum(y[0])-1+lam));
        x0 = sum(y[0]);
        x5 = sum(y[5]);
        for t in range(0,T):
            derivative_2neighboorChoice(y,lam,1,N, derivative);
            y  = y  + h * derivative;
        print('T=',T*conv,'x0=',x0,'x5=',x5);
        conv += 1;
        # for i in range(0,10):
        #     x[t,i] = sum(y[i]);
        #x[t] = sum(y[0]);
    x = array([[0.,0.]]*N);
    for i in range(0,N):
        x[i,0] = i;
        x[i,1] = sum(y[i]);
    savetxt(filename, x);
    return (x[:,1]);

