import math
import copy

def symmetricCheck(matrix):
    '''This function checks if a square matrix is equal to its transpose by comparing elements above and below the diagonal.
    The matrix is considered symmetric if the upper triangle matrix[i][j] equals the lower triangle matrix[j][i] and returns
    a bool 'True' if the matrix is symmetric and 'False' if the matrix is not.'''
    size = len(matrix) #calculates the size of the matrix and stores in 'size'
    for i in range(size): #begins a loop that goes over each row of the matrix and 'i' takes on values 0 to size-1
        for j in range(i + 1, size): #begins another loop that goes over the columns of each matrix and takes on values i+1 to size-1
            if matrix[i][j] != matrix[j][i]: #checks if the matrix is not equal
                return False
    return True

def determinant(matrix):
    '''this function uses a recursive approach to calculate the determinant of a square matrix of any size
    by using a direct formula for a matrix size of 2x2 as the base. Using cofactor expansion for larger
     matrices the determinant is calculated by expanding along the first row'''
    if len(matrix) == 2: #checks if matrix is 2x2 and used as a base for recursion
        return matrix[0][0] * matrix [1][1] - matrix [0][1] * matrix[1][0] #

    det = 0 #det will hold the value of the determine as the function does its iterations
    for column in range(len(matrix)): #begins loop that iterations over each column of the first row matrix
        minor = [row[:column] + row[column+1:] for row in matrix [1:]] #calculates the minor used in the matrix
        det += ((-1) ** column) * matrix[0][column] * determinant(minor)
    #expands the determinant along the first row of matrix through cofactor expansion
    return det #returns the determinant of the matrix

def check_positiveDefinite(matrix):
    '''checks if matrix is positive definite'''
    size = len(matrix) #calculates the size of matrix
    for i in range(1, size + 1): #begins loop that iterates over the values from 1 to size
        submatrix = [row[:i] for row in matrix[:i]] #constructs square submatrix of size 'i*i'
        if determinant(submatrix) <= 0: #calculates the determinant of the submatrix using the 'determinant' function
            return False  #if submatrix is non positive
    return True #if submatrix is positive


def cholesky_decomposition(A):
    '''performs the Cholesky decomposition of a symmetric positive definite matrix A and
    returns the lower triangular matrix L for A = LL^T.'''
    size = len(A) #calculates size of matrix A
    L = [[0.0] * size for _ in range(size)] #creates empty lower triangle matrix L with same dimension as A

    for i in range(size): #begins loop to iterate over rows of A that will fill in the values of L
        for j in range(i + 1): #iterates over columns for lower triangle values
            sum = sum(L[i][k] * L[j][k] for k in range(j)) #sums the products for the corresponding elements of L

            if i == j:  #computes the square root for diagonal elements on the diagonal
                L[i][j] = (A[i][i] - sum) ** 0.5
            else:
                L[i][j] = (A[i][j] - sum) / L[j][j] #computes off diagonal elements of L
    return L


def forward_substitution(L, b):
    '''function to solve the lower triangular system Ly = b for y'''
    size = len(b) #calculates size of vector b
    y = [0 for _ in range(size)] #creates empty vector y for Ly = b

    for i in range(size):
        y[i] = (b[i] - sum(L[i][j] * y[j] for j in range(i))) / L[i][i] #solves using forward substitution

    return y


def backward_substitution(LT, y):
    '''function to solve the upper triangular system L^Tx = y for x accepting non transposed U if necessary'''
    size = len(y)
    x = [0] * size
    for i in range(size - 1, -1, -1):
        x[i] = (y[i] - sum(LT[i][j] * x[j] for j in range(i + 1, size))) / LT[i][i] #solves using backwards substitution
    return x


def cholesky_solve(A, b):
    '''function to solve Ax = b using the Cholesky decomposition'''
    L = cholesky_decomposition(A) #retrives lower triangular matrix L
    LT = [list(row) for row in zip(*L)]  # Transpose L to get LT
    y = forward_substitution(L, b) #solves Ly = b
    x = backward_substitution(LT, y) #solves L^Tx = y
    return x


def LUFactorization(A): #copied and modified from SP2024_Week2
    n = len(A)
    U = [[0.0 for _ in range(n)] for _ in range(n)]
    L = [[0.0 if i != j else 1 for j in range(n)] for i in range(n)]

    for j in range(n):
        for i in range(j+1):
            sum_u = sum(L[i][k] * U[k][j] for k in range(i))
            U[i][j] = A[i][j] - sum_u

        for i in range(j, n):
            sum_l = sum(L[i][k] * U[k][j] for k in range(j))
            if U[j][j] == 0:  # Check to prevent division by zero
                print("Division by zero detected. Matrix might be singular or require pivoting.")
                return L, U  # Early return; Handle this case appropriately.
            L[i][j] = (A[i][j] - sum_l) / U[j][j]

    return L, U

def main_solve(A, b):
    if symmetricCheck(A) and check_positiveDefinite(A): #checks if A is symmetric and positive definite
        print("Using Cholesky method.")
        x = cholesky_solve(A, b)
    else:
        print("Using Doolittle method.")

        L, U = LUFactorization(A) #performs LU Factorization
        y = forward_substitution(L, b) #solves Ly = b for y
        x = backward_substitution(U, y)  #solves Ux = y for x

    print("Solution vector x:", x)
    return x

if __name__ == "__main__":
    #equation set 1
    A1 = [[1, -1, 3, 2], [-1, 5, -5, -2], [3, -5, 19, 3], [2, -2, 3, 21]]
    b1 = [15, -35, 94, 1]
    x1 = main_solve(A1, b1)
    print("Solution for equation set 1 is:", x1)

    #equation set 2
    A2 = [[4, 2, 4, 0], [2, 2, 3, 2], [4, 3, 6, 3], [0, 2, 3, 9]]
    b2 = [20, 36, 60, 122]
    x2 = main_solve(A2, b2)
    print("Solution for equation set 2 is:", x2)