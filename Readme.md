# NTRU Post-Quantum Cryptography Implementation
This project implements the NTRU(N-th degree Truncated polynomial Ring Units) PQC encryption scheme in Python. It allows for the generation of public/private key pairs, encryption of text messages into polynomial ciphertexts, and decryption back to plaintext using lattice-based mathematical principles.

## Introduction
* The ring of polynomials in x with coefficients in $\mathbb{Z}$ is

$$ \mathbb{Z}[x] = \[ \sum_{i=0}^{ \infty } a_i x^i : a_i \in \mathbb{Z},  a_i = 0 \text{ for all but finity many} \] $$

* $\mathbb{Z}_q$ is the finite field where $q$ is prime or power of prime.

$$ \mathbb{Z}_q[x] = \[ \sum_{i=0}^{ \infty } a_i x^i : a_i \in ℤ_q,  a_i = 0 \text{ for all but finity many} \] $$

* N-th degree Truncated Polynomial Ring (of rank N)

$$R = \mathbb{Z}[x]/(x^N -1) = \[ a_0 + a_1 x + \cdots + a_{N−1}x^{N−1} : a_0 , a_1 , … , a_{N−1} \in \mathbb{Z} \] $$

* N-th Degree Truncated Polynomial Ring over a finite field $ℤ_q$

$$R_q = \mathbb{Z}_q[x]/(x^N -1) = \[ a_0 + a_1 x + \cdots + a_{N−1}x^{N−1} : a_0 , a_1 , … , a_{N−1} \in \mathbb{Z}_q \] $$

#### Definition : Center-lift
For any polynomial $f(x) \in R_q$ center-lift of $f$ means shifting coefficients from $[0, q-1]$ to $[-\frac{q}{2} +1 , \frac{q}{2}]$

#### Definition : $T(d_1, d_2)$
$T(d_1, d_2)$ = { $f(x) \in R_q$ : $f(x)$  has  $d_1$  many coefficients equal to 1, $d_2$ many coefficients equal to −1, and the rest are all 0 }.<br>
Polynomials in $T(d_1, d_2)$ are called ternary polynomials.

## Parameters Overview

Parameter | Description
:---:     |  :---
N         | The degree of the polynomial ring. Must be a prime number to ensure the ring $\mathbb{Z}[x]/(x^N-1)$ has a higher probability of producing invertible elements.
p         | The small odd modulus.
q         | The large modulus, Usually a power of 2.
d         | Parameter for the number of 1s and -1s in the ternary polynomials

satisfying $q > (6d + 1)p, \quad gcd(N, q) = 1 = gcd (p, q),\quad 2d + 1 < N$

## Implementation
- `poly_mul(poly1, poly2, N, mod)`: multiplies two polynomials in the ring $R = \mathbb{Z}[x] / (x^N - 1)$. It uses a convolution followed by a reduction step where $x^N$ is treated as $1$, $x^{N+1}$ is treated as $x$ and so on.

- `invert_poly(f, N, p)`: using Extended Euclidean Algorithm (EEA) find the inverse of $f$ modulo $p$, the private key $f_p$.

- `invert_poly_pow2(f, N, q)`: using Hensel Lifting to find the inverse of $f$ modulo 2, then lifting to  $f$ mod $q$, $f_q$.

- `center_lift(poly, mod)`: shifts coefficients from the range $[0, mod-1]$ to $[-\frac{mod}{2}+1, \frac{mod}{2}]$. This "centers" the values, which is the secret to error-free decryption in lattice-based math.

- `generate_ternary_poly(N, num_ones, num_neg_ones)`: generates a random polynomial with a fixed number of 1s and -1s called ternary polynomial.

### Key Generation:
`generate_keys(N, p, q, d)`:
* choose a random polynomial $f(x) \in T(d+1,d)$ and $g(x) \in T(d,d)$
* compute $f_q(x)= (f(x))^{-1} \pmod q$ and $f_p(x)= (f(x))^{-1} \pmod p$ (If either inverse fails to exist, discards this f(x) and choose another one).
* Compute $h(x) = pf_q(x) ∗ g(x) \pmod q$ i.e., $h(x) \in R_q$.

Now $f, f_p$ are the private keys and $h$ is the public key.

### Encryption & Decryption:
* `encrypt(m, h, N, q, d)`:
    * Generates a random "blinding" polynomial $r$.
    * Computes $e = (r \cdot h + m) \pmod q$ where is $m$ is the message polynomial.
* `decrypt(e, f, f_p, N, p, q)`:
    * Calculates $a = (f \cdot e) \pmod q$
    * center-lift $a \pmod p$ to $R$, i.e., coefficients lies in $[-\frac{q}{2}+1, \frac{q}{2}]$.
    * Compute $c = (f_p \cdot b) \pmod p$.
    * center-lift $c \pmod p$ to $R$, i.e., coefficients lies in $[-\frac{p}{2}+1, \frac{p}{2}]$.

### Data Handling:
- `text_to_poly(text, N)`: Converts a string of characters into a binary bitstream, mapping each bit to a polynomial coefficient.

- `poly_to_text(poly, N)`: Reconstructs the 8-bit character chunks from the polynomial coefficients.

# Input & Output
### Input Format:
Input is read from `input.txt` file.
* First line: integer `N`(prime number).
* Second line: integer `p`(small odd modulus).
* Third line: integer `q`(large modulus usually power of 2).
* Forth line: integer `d`(parameter for the number of 1s and -1s in the ternary polynomials).
* Fifth line: Message string(in English).

### Output Format:
If $N$ is small prime ($<30$) then print
- Original Message(in English)
- Message Polynomial
- private key $f$ and $f_p$
- public key $h$
- Encrypted Polynomial
- Random Binding Polynomial
- Decrypted Polynomial
- Recovered Text(in English)

Also print the Success/Failure message.

# How to Run
#### Prerequisites
- Python 3.x
- NumPy (pip install numpy)

### For macOS/Linux
Navigate to the project directory, then run the following commands
```
> source venv/bin/activate
> python3 ntru.py
```
`source venv/bin/activate` for activating the virtual environment, then run using `python3 ntru.py`.

### For Windows
Navigate to the project directory, then run the following commands
```
> venv\Scripts\activate
> python3 ntru.py
```
`venv\Scripts\activate` for activating the virtual environment, then run using `python3 ntru.py`.
