from numpy import sqrt, array, empty
#import matplotlib.pyplot as plt

## Parameters
G = 6.67e-11
M = 1.9891e30 # Mass of the sun (kg)
t0 = 0.0 # intial time (s)
tf = 2e9 # final time (s)
H = tf # Initial step size, which is the ENTIRE interval (s)
n_max = 8 # Maxium allowed no. of splits of the current (sub)interval of
          # integration
delta =  1000 / (3600*24*365.25) # accuracy in position per unit of time (m/s)

# Initial conditions, i.e position of comet in the xy-plane
x0 = 4e12 # m
y0 = 0.0
vx0 = 0.0
vy0 = 500.0  # (m/s)
r0 = array([x0, vx0, y0, vy0], float) # Initial state vector

## Returns the components of the current acceleration and velocity vectors
def f(r):
    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]
    dist = sqrt(x ** 2 + y ** 2)
    ax = -G * M * x / dist ** 3
    ay = -G * M * y / dist ** 3
    return array([vx, ax, vy, \
                     ay], float)


#Lists for storing positions at different time t
xpoints = [x0]
ypoints = [y0]
tpoints = [t0]
# =============================================================================
# The step function for integration of the comet's orbit. It should call itself
# recursively with time step H/2 iff the ticker n > n_max, i.e. if the current
# timestep H needs to be divided into more than 8 parts.
# =============================================================================
def step(r,t,H):
    # Starting modified midpoint step of size H
    n = 1 # We are only taking ONE step initially in the modified
          # midpoint rule
    r1 = r + 0.5*H*f(r)
    r2 = r + H*f(r)

    # Estimate the first row in Richardson extrapolation, with error O(H^2)
    Rnew = empty([1,4],float)
    Rnew[0] = 0.5*(r1 + r2 + 0.5*H*f(r2))

# =============================================================================
# This is where the recursion starts. Rnew will contain the positions at
# t + H, and these are the only points will be stored in memory!
# =============================================================================
    error = 2*H*delta
    while error > H*delta: # H*delta is the tolerance in error for stepsize H
        if n < n_max:
            n += 1
            h = H/n # We're now in the next row of the Richardson extrapolation
                    # Shecme(e.g row 2 for first time around)

            # Starting modified midpoint step of size h_n = H/n
            r1 = r + 0.5*h*f(r)
            r2 = r + h*f(r)
            for i in range(n - 1):
                r1 += h*f(r2)
                r2 += h*f(r1)
# =============================================================================
# Richardson extrapolation starts here in earnest. I've checked that it works
# for easier dynamical systems that do NOT demand any recursion.
# =============================================================================

            Rold  = Rnew # Making room for a redefinition of the current row
            Rnew  = empty([n,4], float)
            Rnew[0] = 0.5*(r1 + r2 + 0.5 * h * f(r2))
            for m in range(1,n):
                # Error estimate to O(h^(2m))
                epsilon = (Rnew[m-1] - Rold[m-1]) / ((n / (n - 1))**(2*m) - 1)
                Rnew[m] = Rnew[m-1] + epsilon
            # "Euclidian error" estimate
            error = sqrt(epsilon[0] ** 2 + epsilon[2] ** 2)
        else: break
# =============================================================================
# THIS is the part that drives me nuts! If we still haven't reached our goal,
# Then subdivide the current interval [t, t+H) into two parts,
# [t, H/2), and [T + H, T + H/2). How to append lists and return positions
#in the correct order??
# =============================================================================
    if error > H * delta and n > n_max:
        # What to do with these recursive calls??
        R1 = step(r, t, H / 2)
        R2 = step(r, t + H/2, H/2)
    else:
        r += Rnew[n-1]
        xpoints.append(Rnew[n-1, 0]) # Best estimate of x at t+h
        ypoints.append(Rnew[n-1, 2]) # --------||------ y at t+h
        tpoints.append(t+h)

# =============================================================================
# End of my misery...
# =============================================================================

# calling the functions with initial conditions and the "one-step" stepsize
step(r0, t0, H)


