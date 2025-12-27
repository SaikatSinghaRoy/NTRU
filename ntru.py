# N – prime number, the polynomials in the truncated polynomial ring have degrees atmost N-1
# q – large modulus
# p – small odd modulus
# d - Parameter for the number of 1s and -1s in the ternary polynomials
# satisfying  q > (6d + 1)p, gcd(N, q) = 1 = gcd (p, q), 2d + 1 < N


import numpy as np

########################
### Helper Functions ###
########################
def poly_add(f, g, mod):
    return np.mod(f + g, mod)

def poly_mul(f, g, N, mod):
    # Standard polynomial multiplication
    res = np.convolve(f, g)
    # Reduce by x^N - 1 (x^N = 1, x^(N+1) = x,  and so on)
    final_res = np.zeros(N)
    for i, val in enumerate(res):
        final_res[i % N] += val
    return np.mod(final_res, mod).astype(int)

def poly_division_mod(f, g, p):

    f = np.trim_zeros(f, 'b')
    g = np.trim_zeros(g, 'b')
    
    deg1 = len(f) - 1
    deg2 = len(g) - 1
    
    if deg2 < 0: raise ZeroDivisionError()
    if deg1 < deg2: return np.array([0]), f
    
    rem = np.copy(f) % p
    quot = np.zeros(deg1 - deg2 + 1, dtype=int)
    
    # inverse of leading coefficient of divisor poly
    inv_lead = pow(int(g[-1]), -1, p)
    
    for i in range(deg1 - deg2, -1, -1):
        if len(rem) < i + deg2 + 1: continue
        factor = (rem[i + deg2] * inv_lead) % p
        quot[i] = factor
        for j in range(deg2 + 1):
            rem[i + j] = (rem[i + j] - factor * g[j]) % p
            
    return quot % p, np.trim_zeros(rem, 'b') % p

def inv_of_poly(f, N, p):

    # Define x^N - 1
    modulus = np.zeros(N + 1, dtype=int)
    modulus[0] = -1
    modulus[N] = 1   # modulus = x^N -1
    
    # old_r = f * old_s + modulus * something
    old_r, r = (f % p), (modulus % p)
    old_s, s = np.array([1]), np.array([0])
    
    while np.any(r):
        quot, rem = poly_division_mod(old_r, r, p)
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

def inv_of_poly_q(f, N, q):
    # Find inverse mod 2 using the EEA function
    try:
        f_inv = inv_of_poly(f, N, 2)
    except ValueError:
        raise ValueError("f is not invertible mod 2, so it cannot be lifted to mod q")

    # Lift mod 2 to mod 4, 8, 16... up to q
    current_mod = 2
    while current_mod < q:
        current_mod *= current_mod
        # f_inv = f_inv * (2 - f * f_inv) mod current_mod
        
        #  (f * f_inv) 
        f_prod = poly_mul(f, f_inv, N, current_mod)
        
        # Calculate (2 - f_prod)
        poly_two = np.zeros(N, dtype=int)
        poly_two[0] = 2
        diff = np.mod(poly_two - f_prod, current_mod)
        
        f_inv = poly_mul(f_inv, diff, N, current_mod)
        
    return np.mod(f_inv, q).astype(int)

def generate_ternary_poly(N, num_ones, num_neg_ones):
    x = np.zeros(N, dtype=int)
    indices = np.random.choice(N, num_ones + num_neg_ones, replace=False)
    x[indices[:num_ones]] = 1
    x[indices[num_ones:]] = -1
    return x

def center_lift(x, q):
    # change coefficients from [0, q-1] to [-q/2+1, q/2] 
    new_x = np.copy(x)
    for i in range(len(new_x)):
        val = new_x[i] % q
        if val > q // 2:
            val -= q
        new_x[i] = val
    return new_x

#######################
###  Data Handling  ###
#######################
def text_to_poly(text, N):
    
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
    
    # Convert coefficients back to a bit string
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
    
    while True:
        try:
            f = generate_ternary_poly(N, d+1, d)
            g = generate_ternary_poly(N, d, d)

            # Compute f_p (inverse of f mod p)
            f_p = inv_of_poly(f, N, p)
            # Compute f_q (inverse of f mod q)
            f_q = inv_of_poly_q(f, N, q)
            break

        except ValueError:
            continue # when f in not invertible try new f

    # Public Key h = (p * f_q * g) mod q
    pf_q = np.mod(p * f_q, q)
    h = poly_mul(pf_q, g, N, q)
    
    return (f, f_p), h

### Encryption
def encrypt(m, h, N, q, d):
    
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

    # a = (f * e) mod q
    a = poly_mul(f, e, N, q)
    
    # Center-lift a modulo q
    a_lifted = center_lift(a, q)
    
    # c = (f_p * a_lifted) mod p
    c = poly_mul(f_p, a_lifted, N, p)
    
    # center lift c (mod p)
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
    print(f"\nOriginal Message:        {poly_to_text(m, N)}")
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
