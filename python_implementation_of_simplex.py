import numpy as np 
import math


"""Theory:
type:
max z=cx
Ax<=b
x>=0

z = cb*xb + cn*xn
Splitting A also into basis matrix B and Non-basis matrix N we will have
B*xb + N*xn = b
x>=0

now B*xb = b - N*xn
xb = B_inv*b - B_inv*N*xn

z = cb*xb + cn*xn
replacing xb we have

z = cb*(B_inv*b) - (cb*B_inv*N - cn)*xn

Step 0:
Determine initial BFS

Optimality test:
Evaluate Zj - Cj for each non-basic varialble

Zj-Cj = cb*B_inv* Nj - cj

if all (Zj-Cj) > 0     ---> Solution is optimal
if min(Zj-Cj) = 0      --->Multiple optimal solution
if min(Zj-Cj) < 0      --->Iterate further

Choosing leaving basic variable:
if leaving_bv is found successfully :
    continue iterating
else:
    Problem is unbounded


"""


def simplex(c,A,b):
    """This function runs iteration for revised simplex method
        
       Problem :
       Maximize cx\n
       subject to Ax = b\n
       x>=0

       """
    # converting all into numpy arrays if such is not the case
    c = np.array(c)
    A = np.array(A)
    b = np.array(b)


    print("Coefficients of objective function ")
    print(c)
    print("The constraint matrix with slacks is")
    print(A)
    print("The RHS of constraints is")
    print(b,end = '\n\n\n')


    # initializing B
    xb,B,xn,N = initialize_B(A)
    B_inv = np.linalg.inv(B)
    
    z,RHS = Value_finder(B,b,c,xb)

    print("Initial basic solution:")
    print("The basic variables are {}".format(xb))
    print("The values of basic variables are\n {}".format(RHS))
    print("THe value of objective function is {}".format(float(z)))

    optimality,entering_nb = optimality_test(c,xb,xn,A,B)


    while optimality == False:
        
        print("Optimality Test = {}".format(optimality))
        print("The entering Non-Basic variable is {}".format(entering_nb))

        leaving_bvar = min_ratio_test(entering_nb,xb,A,RHS,B_inv)

        # if minimum_ratio_test fails the problem is unbounded
        if leaving_bvar == -1:
            print("No leaving basic variable")
            print("The problem is unbounded")
            return (-1000,[],[])    # randomly setted flag for unbounded problem
        else:
            print("The leaving basic variable is {}".format(leaving_bvar))
        
        # else change the basis
        xb.remove(leaving_bvar)
        xn.append(leaving_bvar)
       
        xb.append(entering_nb)
        xn.remove(entering_nb)

        # get the matrices B and N
        B = A[ : ,xb]
        N = A[ : ,xn]

        # evaluate B_inv
        B_inv = np.linalg.inv(B)

        # find objective value and values for basic variables        
        z,RHS = Value_finder(B,b,c,xb)

        print("\n\nAfter new iteration :")
        print("Basic variables = \n{}".format(xb))
        print("Value of basic variables =\n {}".format(RHS))
        print("Value of Objective function z = {}".format(float(z)))

        # check for optimality
        optimality,entering_nb = optimality_test(c,xb,xn,A,B)

    if entering_nb == -2:
        print("Multiple optimal solutions as one non-basic variable has now got 0 objective coefficient")
        

    # when optimal solution has been found
    print("\n\nOptimality has been reached")
    print("The basic variables are")
    print(xb)
    print("Value of basic variables =")
    print(RHS)
    print("The optimal value of objective function is {}".format(float(z)))

    return (z,xb,RHS)

    


# tested - working correctly
def initialize_B(A):
    """Takes in matrix A and returns matrix B for initial BFS"""
    # note that we make 0 to be initial BFS

    dims = np.shape(A)

    m = dims[0] # number of constraints
    n = dims[1] # number of variables including slack


    xb = []
    xn = []

    for index in range(0,n):
        if index < (n-m):
            xn.append(index)
        else:
            xb.append(index)      


    B = A[:,xb]
    N = A[:,xn]
    
    return (xb,B,xn,N) 

# tested - working correctly
def optimality_test(c,xb,xn,A,B):
    """Returns true if solution is optimal
       returns the greatest cj-zj otherwise"""
    
    # getting coefficients of Basic variables
    cb = c[xb] # note that c must be numpy array
    B_inv = np.linalg.inv(B)

    # evaluating matrix multplication of cb * B_inv
    cb_binv_prd = np.matmul(cb,B_inv)


    zj_min = 1          # minimum of all zj-cj, initialized to positive value
    entering_nb = -1    #initialize with invalid value

    # iterate through non-basic variables
    for index in xn:
        # getting column of non-basic variable        
        Nj = A[:,index]
        # evaluating zj-cj
        potential = np.matmul(cb_binv_prd,Nj) - c[index]

        if potential < zj_min:
            # update entering_nb
            entering_nb = index
            zj_min = potential
        
    # retruns the index of entering NBV if solution is not optimal
    # otherwise returns -1 to show optimal solution found

    if entering_nb == -1:
        return True,-1
    elif (zj_min==0):
        return True,-2
    else:
        return False,entering_nb

# tested - working correctly
def min_ratio_test(entering_nb,xb,A,RHS,B_inv):
    """We need to minimum positive (B_inv*b)jth/ (B_inv*Nj)"""

    leaving_index = -1 # initialize to invalid value

    # extracting column of entering Non-Basic Variable
    Nj = A[:,entering_nb:entering_nb+1]     # note this is done to prevent dims like (2,?)

    
    enterinng_column = np.matmul(B_inv, Nj )
    
    mrt_array = RHS/enterinng_column

    # getting the minimum ratio
    mrt_min = np.min(mrt_array)

    # if minimimum ratio is positive update values
    if mrt_min > 0:
        mrt_index = np.argmin(mrt_array)
        # finding true index in original A
        leaving_index = xb[mrt_index]

    # if MRT fails return -1
    # else return leaving_basic variable
    return leaving_index
           

# tested - working correctly
def Value_finder(B,b,c,xb):
    """Returns objective function value z along with values of basic variables
       inputs: B,b,c,xb(indices of basic variables)
       ouptus: z, RHS (values of basic variabless"""
    
    B_inv = np.linalg.inv(B)
    
    cb = c[xb]  # cefficients of basic variables
    RHS = np.matmul(B_inv,b) #values of basic variables

    z = np.dot(cb,RHS)

    return z,RHS

# test run case below
def main():
    # testing starts here

    A = [[30,20,1,0],
         [5,10,0,1]]

    b = [ [300],
          [110] ]

    c = [6,8,0,0]

    z,xb,RHS=simplex(c,A,b)
    


    return 0

main()
        

