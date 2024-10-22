# Binomial-Asset-Pricing-Model-
A Python implementation of various binomial asset pricing models, including CRR, JR, EQP, and TRG, with a comparison to the Black-Scholes model for option pricing.

We treat the binomial tree as a network with nodes `[i,j]` with `i` representing the time steps and `j` representing the number of ordered prive outcome. Lowest to Highest - starts from bottom of the tree.

The following methods are implemented within the `binomial_tree_slow` algorithm:

- **Cox, Ross, and Rubinstein (CRR) Method**: Chooses equal jump sizes.
- **Jarrow and Rudd (JR) Method**: Chooses equal risk-neutral probabilities.
- **Equal Probabilities (EQP) Method**: Uses equal risk-neutral probabilities under a logarithmic asset pricing tree.
- **Trigeorgis (TRG) Method**: Chooses equal jump sizes under a logarithmic asset pricing tree.

These methods also work for the `binomial_tree_fast` algorithm.
***

## Installation of Dependencies

To run the code, you need to install the following Python libraries:

```pip install numpy scipy matplotlib py_vollib```
***
# Binomial Tree Representation
Stock tree can be represented using nodes `[i,j]` and the initial stock price $S_0$
  $S_{i,j}$ $=$ $S_0u^jd^{i - j}$
  $C_{i,j}$ represents contract price at each node `[i,j]`, where $C_{N,j}$ represents final payoiff function that we can define.

*For this project, we will price an European Call, so* $C_{N,j}$ $=$ $max(S_{N,j} - K,0)$

```python
S0 = 100      # Initial stock price
K = 110       # Strike price
T = 0.5       # Time to maturity in years
r = 0.06      # Annual risk-free rate
N = 100       # Number of time steps
sigma = 0.3   # Annualized stock price volatility
opttype = 'C' # Option type 'C' for Call, 'P' for Put
```
***
## Methods Overview
- Cox, Ross, and Rubinstein (CRR) Method
This method assumes equal jump sizes and calculates option prices using the following formula:
```python
def CRR_method(K, T, S0, r, N, sigma, opttype='C'):
    # Precompute constants
    dt = T/N
    u = np.exp(sigma * np.sqrt(dt))
    d = 1/u
    q = (np.exp(r * dt) - d) / (u - d)
    disc = np.exp(-r * dt)

    # Initialize asset prices at maturity
    S = np.zeros(N+1)
    S[0] = S0 * d**N
    for j in range(1, N+1):
        S[j] = S[j-1] * u/d

    # Initialize option values at maturity
    C = np.zeros(N+1)
    for j in range(0, N+1):
        C[j] = max(0, S[j] - K) if opttype == 'C' else max(0, K - S[j])

    # Step backward through tree
    for i in np.arange(N, 0, -1):
        for j in range(0, i):
            C[j] = disc * (q * C[j+1] + (1-q) * C[j])

    return C[0]
```
- Similar methods are provided for:

1. Jarrow and Rudd (JR) Method
2. Equal Probabilities (EQP) Method
3. Trigeorgis (TRG) Method


## Comparison of Methods
The performance of these methods can be compared by evaluating the convergence of the option price as a function of the number of time steps (N). We use the Black-Scholes formula for reference and compare the following methods:

- Cox, Ross, and Rubinstein (CRR)
- Jarrow and Rudd (JR)
- Equal Probabilities (EQP)
- Trigeorgis (TRG)

```python
import matplotlib.pyplot as plt
from py_vollib.black_scholes import black_scholes as bs

CRR, JR, EQP, TRG = [], [], [], []

periods = range(10, 500, 10)
for N in periods:
    CRR.append(CRR_method(K, T, S0, r, N, sigma, opttype='C'))
    JR.append(JR_method(K, T, S0, r, N, sigma, opttype='C'))
    EQP.append(EQP_method(K, T, S0, r, N, sigma, opttype='C'))
    TRG.append(TRG_method(K, T, S0, r, N, sigma, opttype='C'))

BS = [bs('c', S0, K, T, r, sigma) for _ in periods]

plt.plot(periods, CRR, label='Cox_Ross_Rubinstein')
plt.plot(periods, JR, label='Jarrow_Rudd')
plt.plot(periods, EQP, label='EQP')
plt.plot(periods, TRG, 'r--', label='Trigeorgis')
plt.plot(periods, BS, 'k', label='Black-Scholes')
plt.legend(loc='upper right')
plt.show()

```
***
## Conclusion

This project showcases the implementation and comparison of several binomial asset pricing models (CRR, JR, EQP, and TRG) against the Black-Scholes model. Key insights include:

1. **Convergence**: All models converge to the Black-Scholes price as the number of time steps (N) increases, with CRR and JR converging the fastest.
2. **Performance**: CRR and JR are more accurate and computationally efficient, especially for small N, while EQP and TRG require more time steps to converge due to their complexity.
3. **Applicability**: CRR and JR are ideal for quick, accurate pricing, while EQP and TRG offer more flexibility but at a higher computational cost.

This comparison helps users choose the right model based on their accuracy and computational requirements.
