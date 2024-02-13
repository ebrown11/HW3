import math
def gamma_Function(z):
    '''here is the Lanczos approximation to compute the gamma function numerically'''
    coef = [1.000000000190015, 76.18009172947146, -86.50532032941677, 24.01409824083091,
            -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5] #coefficients for Lanczos approx.
    g = 7 #constant for lanczos for accuracy
    z -= 1 # formula requires z-1
    x = 1.000000000190015 #starting value of approximation
    for i in range(len(coef)): #begins loop over the coef. starting from the second coef.
        x += coef[i] / (z + i + 1) #calculates the sum portion
    t = z + g + 0.5 #used for final exponentiation and multiplication
    return (2.5066282746310007 * x) * (t ** (z + 0.5)) * (math.exp(-t))

def K_Function(m):
    '''calculates the constant K_m for the t-distribution'''
    return (gamma_Function(.5 * m + .5)) / (math.sqrt(m * math.pi) * gamma_Function(.5 * m))

def t_Distribution(u, m):
    '''calculates F(z) for a given t-value and degrees of freedom m using the numerical in '''
    return (1 + u**2 / m) ** (-(m + 2) / 2)

def simpson_Rule(a, b, f, n):
    '''function that performs numerical integration of the function f over the interval using Simpson's rule and
    returns the approximate value of the integral of f from a to b'''
    if n % 2 != 0: #checks if n is not even
        raise ValueError("n must be an even number.")

    h = (b - a) / n #calculates width of each sub interval
    total = f(a) + f(b)

    for i in range(1, n): #iterates from 1 to n-1
        x = a + i * h
        if i % 2 == 0:
            total += 2 * f(x)
        else:
            total += 4 * f(x)

    return (h / 3) * total
def Fz(x, m, n):
    '''calculates F(s) up to x for the t-distribution with m degrees'''
    f = lambda u: t_Distribution(u, m)

    lower_limit = -100

    result = simpson_Rule(lower_limit, x, f, n)
    return result

def user_Input():
    '''function that asked for z value and degrees of freedom and storing as a float and integer, respectively'''
    while True: #begins loop until valid inputs are accepted
        try:
            degrees = int(input("Enter the degrees of freedom as a whole number: "))
            if degrees <= 0:
                print("Degrees of freedom must be positive. Try Again")
                continue

            z_value = float(input("Enter the z-value: "))
            break
        except ValueError:
            print("Invalid input.")
    return degrees, z_value

def main():
    n = 1000 #number of intervals for simpson's rule
    degrees, z_value = user_Input() #calls user_Input function for freedom of degrees and z value
    result = Fz(z_value, degrees, n) #calculates F(z)
    print(f"The cumulative probability F(z) up to {z_value} with {degrees} degrees of freedom is: {result}")

if __name__ == "__main__":
    main()
