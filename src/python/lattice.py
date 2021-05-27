import numpy as np

def update(x):
    raise NotImplementedError
    for j in range(0, N):
        old_x = x[j]
        old_Sj = S(j, x)
        x[j] = x[j] + uniform(-eps, eps)
        dS = S(j, x) - old_Sj
        if dS > 0 and exp(-dS) < uniform(0, 1):
            x[j] = old_x

def S(j, x):
    jp = (j + 1) % N
    jm = (j - 1) % N
    raise NotImplementedError
    return (a/2) * x[j]**2 + (x[j]/a) * (x[j] - x[jp] - x[jm])

def compute_G(x, n):
    raise NotImplementedError

def monte_carlo_average(x, G):
    raise NotImplementedError

def bootstrap(G):
    raise NotImplementedError

def bin(G, bin_size):
    raise NotImplementedError

def avg(G):
    raise NotImplementedError

def sdev(G):
    raise NotImplementedError

def deltaE(G):
    raise NotImplementedError

def bootstrap_deltaE(G, nbstrap=100):
    raise NotImplementedError

def main():
    raise NotImplementedError


if __name__ == "__main__":
    main()
