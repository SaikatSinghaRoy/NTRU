# N – prime number, the polynomials in the truncated polynomial ring have degrees atmost N-1
# q – large modulus
# p – small odd modulus
# d - Parameter for the number of 1s and -1s in the ternary polynomials
# satisfying  q > (6d + 1)p, gcd(N, q) = 1 = gcd (p, q), 2d + 1 < N


import numpy as np

########################
### Helper Functions ###
########################
def poly_add(poly1, poly2, mod):
    return np.mod(poly1 + poly2, mod)

def poly_mul(poly1, poly2, N, mod):
    """Multiplies two polynomials in the ring R = Z[x] / (x^N - 1)"""
    # Standard polynomial multiplication
    res = np.convolve(poly1, poly2)
    # Reduce by x^N - 1 (x^N = 1, x^(N+1) = x,  and so on)
    final_res = np.zeros(N)
    for i, val in enumerate(res):
        final_res[i % N] += val
    return np.mod(final_res, mod).astype(int)

def poly_div_mod(poly1, poly2, p):
    """Robust polynomial long division modulo p"""
    poly1 = np.trim_zeros(poly1, 'b')
    poly2 = np.trim_zeros(poly2, 'b')
    
    deg1 = len(poly1) - 1
    deg2 = len(poly2) - 1
    
    if deg2 < 0: raise ZeroDivisionError()
    if deg1 < deg2: return np.array([0]), poly1
    
    rem = np.copy(poly1) % p
    quot = np.zeros(deg1 - deg2 + 1, dtype=int)
    
    # inverse of leading coefficient of divisor poly
    inv_lead = pow(int(poly2[-1]), -1, p)
    
    for i in range(deg1 - deg2, -1, -1):
        if len(rem) < i + deg2 + 1: continue
        factor = (rem[i + deg2] * inv_lead) % p
        quot[i] = factor
        for j in range(deg2 + 1):
            rem[i + j] = (rem[i + j] - factor * poly2[j]) % p
            
    return quot % p, np.trim_zeros(rem, 'b') % p

def invert_poly(f, N, p):
    """ Extended Euclidean Algorithm for Polynomials in Zp[x]/(x^N - 1) """
    # Define x^N - 1
    modulus = np.zeros(N + 1, dtype=int)
    modulus[0] = -1
    modulus[N] = 1   # modulus = x^N -1
    
    # Standard EEA: old_r = f * old_s + modulus * something
    old_r, r = (f % p), (modulus % p)
    old_s, s = np.array([1]), np.array([0])
    
    while np.any(r):
        quot, rem = poly_div_mod(old_r, r, p)
        old_r, r = r, rem
        
        # s_new = old_s - quot * s
        prod = np.convolve(quot, s) % p
        
        # Subtracting polynomials of different lengths
        max_len = max(len(old_s), len(prod))
        new_s = np.zeros(max_len, dtype=int)
        new_s[:len(old_s)] += old_s
        new_s[:len(prod)] -= prod
        old_s, s = s, new_s % p

    # the GCD (old_r) must be a constant (degree 0)
    if len(old_r) > 1:
        raise ValueError("f and x^N-1 are not coprime. Try a different f.")
        
    # Scale old_s by the inverse of the constant GCD
    inv_gcd = pow(int(old_r[0]), -1, p)
    # Reduce result mod x^N - 1
    final_res = np.zeros(N, dtype=int)
    for i in range(len(old_s)):
        final_res[i % N] = (final_res[i % N] + old_s[i] * inv_gcd) % p
    return final_res

def invert_poly_pow2(f, N, q):
    """
    Computes inverse of f mod (x^N - 1) and mod q (where q is power of 2) using Hensel Lifting.
    """
    # Find inverse mod 2 using the EEA function
    try:
        f_inv = invert_poly(f, N, 2)
    except ValueError:
        raise ValueError("f is not invertible mod 2, so it cannot be lifted to mod q")

    # Lift mod 2 to mod 4, 8, 16... up to q
    current_mod = 2
    while current_mod < q:
        current_mod *= current_mod
        # Formula: f_inv = f_inv * (2 - f * f_inv) mod current_mod
        
        # Calculate (f * f_inv) in the ring
        f_prod = poly_mul(f, f_inv, N, current_mod)
        
        # Calculate (2 - f_prod)
        # Note: 2 is the polynomial [2, 0, 0, ...]
        poly_two = np.zeros(N, dtype=int)
        poly_two[0] = 2
        diff = np.mod(poly_two - f_prod, current_mod)
        
        # Update f_inv
        f_inv = poly_mul(f_inv, diff, N, current_mod)
        
    return np.mod(f_inv, q).astype(int)

def generate_ternary_poly(N, num_ones, num_neg_ones):
    """ Generates a random polynomial with a fixed number of 1s and -1s called ternary polynomial """
    poly = np.zeros(N, dtype=int)
    indices = np.random.choice(N, num_ones + num_neg_ones, replace=False)
    poly[indices[:num_ones]] = 1
    poly[indices[num_ones:]] = -1
    return poly

def center_lift(poly, q):
    """ Shifts coefficients from [0, q-1] to [-q/2+1, q/2] """
    new_poly = np.copy(poly)
    for i in range(len(new_poly)):
        val = new_poly[i] % q
        if val > q // 2:
            val -= q
        new_poly[i] = val
    return new_poly

#######################
###  Data Handling  ###
#######################
def text_to_poly(text, N):
    """
    Converts a string to a polynomial of degree N-1.
    Each character is converted to bits, and each bit becomes a coefficient.
    """
    # Convert text to a bit string
    bits = ''.join(format(ord(i), '08b') for i in text)
    
    if len(bits) > N:
        print(f"Warning: Text too long for N={N}. Truncating.")
        bits = bits[:N]
    
    # Initialize polynomial with zeros
    poly = np.zeros(N, dtype=int)
    
    # Map bits to coefficients {0, 1}
    for i, bit in enumerate(bits):
        poly[i] = int(bit)
        
    return poly

def poly_to_text(poly, N):
    """
    Converts a polynomial back into a string.
    Expects coefficients to be 0 or 1.
    """
    # Convert coefficients back to a bit string
    # We use center_lift(poly, 3) first to ensure we have 0s and 1s
    bits = "".join(str(abs(int(c))) for c in poly)
    
    # Group into 8-bit chunks and convert to characters
    chars = []
    for i in range(0, len(bits), 8):
        byte = bits[i:i+8]
        if len(byte) < 8: break # Ignore trailing incomplete bytes
        chars.append(chr(int(byte, 2)))
        
    return "".join(chars)


#######################
### NTRU Core Logic ###
#######################
### Key Generation
def generate_keys(N, p, q, d):
    """
    Generating Keys f, f_p, h :: private keys f and f_p :: public keys h
    """
    while True:
        try:
            f = generate_ternary_poly(N, d+1, d)
            g = generate_ternary_poly(N, d, d)

            # Compute f_p (inverse of f mod p)
            f_p = invert_poly(f, N, p)
            # Compute f_q (inverse of f mod q)
            f_q = invert_poly_pow2(f, N, q)
            break

        except ValueError:
            continue # when f in not invertible try new f

    # Public Key h = (p * f_q * g) mod q
    pf_q = np.mod(p * f_q, q)
    h = poly_mul(pf_q, g, N, q)
    
    return (f, f_p), h

### Encryption
def encrypt(m, h, N, q, d):
    """
    Encrypts a message polynomial m using public key h
    e = r * h + m (mod q)
    m: polynomial with coefficients in {-1, 0, 1} len(m) = N
    h: public key polynomial
    r: random binding polynomial
    """
    # Generate a random binding polynomial r
    r = generate_ternary_poly(N, d, d)
    
    # Calculate e = (r * h + m) mod q

    # r * h mod (x^N - 1)
    rh = poly_mul(r, h, N, q)
    # (rh + m) mod q
    e = np.mod(rh + m, q)
    
    return e, r

### Decryption
def decrypt(e, f, f_p, N, p, q):
    """ Decrypts ciphertext e using private keys f and f_p """
    # Step 1: a = (f * e) mod q
    a = poly_mul(f, e, N, q)
    
    # Step 2: Center-lift a modulo q
    # This recovers the polynomial p * r * g + f * m without the mod q wrapping
    a_lifted = center_lift(a, q)
    
    # Step 3: c = (f_p * a_lifted) mod p
    c = poly_mul(f_p, a_lifted, N, p)
    
    # Step 4: center lift c (mod p)
    return center_lift(c, p)


if __name__ == "__main__":

    # Parameters for NTRU
    with open('input.txt', 'r') as f:
        N = int(f.readline().strip())
        p = int(f.readline().strip())
        q = int(f.readline().strip())
        d = int(f.readline().strip())

        # original text read from input.txt
        original_text = f.readline().strip()

    # convert origin message to polynomial
    m = text_to_poly(original_text, N)
    print(f"\nOriginal Message:        {original_text}")
    if(N < 30):
        print(f"Message Polynomial:        {m}")

    # key generation
    (f, f_p), public_key = generate_keys(N, p, q, d)
    print("\nKey generation: DONE")
    if(N < 30):
        print(f"f =                        {f}")
        print(f"f_p =                      {f_p}")
        print(f"Public Key:                {public_key}")

    # encryption
    encryp_poly, random_poly = encrypt(m,public_key, N, q, d)
    print("\nEncryption: DONE")
    if(N < 30):
        print(f"Encrypted Polynomial:      {encryp_poly}")
        print(f"Random Binding Polynomial: {random_poly}")

    # decryption
    print("\nDecryption: DONE")
    decrypt_poly = decrypt(encryp_poly, f, f_p, N, p, q)

    # convert decrypted polynomial to readable text format
    recovered_text = poly_to_text(decrypt_poly, N)
    if(N < 30):
        print(f"Decrypted Polynomial:       {decrypt_poly}")
        print(f"Recovered Text:             {recovered_text}")
    
    # indication of success
    if np.array_equal(m, decrypt_poly):
        print("\nSuccess! The decryption is perfect.")
    else:
        print("\nDecryption failed.")
