from numpy import *
from numpy.matlib import randn
from numpy.linalg import inv,det

from random import *
from numpy.random import *


def predict(x, P, A, Q, B, u):
    x = dot(A,x) + dot(B, u)   
    P = dot(A, dot(P, A.T)) + Q
    return(x,P)

def update(x, P, z, H, R):
    
    Hx = dot(H, x)
    IS = R + dot(H, dot(P, H.T))
    PH = dot(P,H.T)
    K= dot(PH,inv(IS))
    x = x + dot(K, (z-Hx)[0])#by miec skalar, a nie wektor 1 elementowy ...
 
    P = P - dot(K, dot(H, P)) 
    return (x,P)


def KalmanInit():
    # intial parameters
    n_iter = 7 #iteration number
    #X - state vector = [teta,omega, alfa]
    #teta - angle
    #omega(n)= [teta(n+1)-teta(n-1)]/2  temporary freq.
    #alfa(n)=(teta(n+1)-teta(n))- (teta(n)-teta(n-1))= teta(n+1)-2teta(n)+teta(n-1) modulation factor
    X=zeros((n_iter,3))
    X_estimate=zeros((n_iter,3))
    X_predict=zeros((n_iter,3))
    teta_list=zeros(n_iter+1)  
    omega_list=zeros(n_iter)  
    alfa_list=zeros(n_iter)  
     
    for n in range(n_iter+1):
        teta_list[n]=n+n*n
     
    #temporary freq. and modulation factor - not important if H= [1 0 0]   
  
    #for n in range(1,n_iter):# dla n = 0 -> omega,alfa=0
    #    omega_list[n]=(teta_list[n+1]-teta_list[n-1])/2
    #    alfa_list[n]=teta_list[n+1]-2*teta_list[n]+teta_list[n-1]
    for n in range(n_iter):
        X[n][0]=teta_list[n]
        #X[n][1]=omega_list[n]
        #X[n][2]=alfa_list[n]
     
    X=array(X)
    A  = array([[1, 1, 0.5], [0, 1, 1], [0, 0, 1]])
    H = array([1, 0, 0])
    B = eye(X.shape[1])
    #print B.shape=(3,3)
    U = zeros((n_iter,X.shape[1])) 
    Q = eye(X.shape[1]) 
    R = eye(X.shape[1])
    P_init  = diag((0.01, 0.01, 0.01))    
    Z=zeros((n_iter,1))
    V = normal(0,1,n_iter)

    for n in range(n_iter):
        Z[n]=dot(H,X[n])+V[n]
         
    return n_iter,X_estimate,X_predict,A,H,Q,B,U,R,P_init,Z 
    
def ccmKalman():

    n_iter,X_estimate,X_predict,A,H,Q,B,U,R,P_estimate,Z=KalmanInit()
    print "observation: ", Z
    
    for n in range(1,n_iter):
        (X_predict[n],P_predict)=predict(X_estimate[n-1],P_estimate, A, Q, B, U[n-1]) 
        (X_estimate[n],P_estimate)=update(X_predict[n],P_predict, Z[n], H, R)
        
        
    print "X_estimate = ", X_estimate.T[0]
      

ccmKalman()