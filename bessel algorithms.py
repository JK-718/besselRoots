import scipy.special

#Define bessel function of order zero and its derivative

def f(x):
    return  scipy.special.jv(0,x)

def fprime(x):
    return -(scipy.special.jv(1,x))
    
#Begin code for bisection algorithm to find roots of function
def bisection(a, b, n, e, d):
    c = (a+b)/2
    maxiter = 0
    while abs((b-a)) > d and abs(f((a+b)/2)) > e and maxiter != n :     
        if f(c) == 0:
            exit      
        elif f(a)*f(c) < 0:
            b = c
            maxiter += 1   
        else:
            a = c
            maxiter += 1
        c = (a+b)/2.0
    if maxiter == n:
       print("reached max number of iterations")     
    return c, maxiter

#Secant algorithm to find roots
def secant(a, b, n, e, d):
    maxiter = 0
    while abs(f(b))>e and abs(b - a)>d and maxiter != n:      
        if f(b) - f(a) == 0:
            return       
        else:
            m = b-f(b)*(b-a)/(f(b)-f(a))
            a = b
            b = m
            maxiter += 1   
    if maxiter == n:
        print("reached max number of iterations")    
    return m, maxiter

#Newton's method to find roots
def newton(x, n, e):
    maxiter = 0
    while  maxiter != n and abs(f(x))>e:  
        x = x - f(x)/fprime(x)
        maxiter += 1
        if maxiter == n:
            print("reached max number of iterations")
    return x, maxiter

#a, b are bounds for where to find roots
#used in bisection and secant
a = [1,4,8]
b = [3,6,10]

#max number of iterations
n = 100

#error values
e = [10E-6, 10E-7, 10E-8,10E-9,10E-10]
d = 10E-16

#initial guesses for Newton's method
x = [3, 6, 10]

#loop over the ranges given
#compare accuracies and computation time
#Bounds are for bisection and secant
#Initial guess is for Newton
for i, j, k in zip(a,b,x):
    print("\n Bounds", i, j, "& Initial guess", k)
    for l in e:
        print("Error level", l)
        print("Bisection", bisection(i, j, n, l, d))
        print("Secant", secant(i, j, n, l, d))
        print("Newton", newton(k, n, l))
        

