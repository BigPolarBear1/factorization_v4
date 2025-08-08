# factorization_v4

For the current fastest version, look at v3 (which is also still a work in progress as I rewrite that project to C).

This is however the cutting edge of my number theoretical research and the direction I should continue in.

The current PoC creates a hashmap that shows for every linear coefficient mod p at which quadratic coefficient mod p it is found.
I.e if we have as valid solution 6x<sup>2</sup>+4x = 0 mod N then we will find 6 as a valid solution as quadratic coefficient for linear coefficient 4 in EVERY MOD P, without exception.
This "idea" is similar to using legendre symbols to determine squaredness of an integers by checking if its a square mod multiple primes and the more primes moduli it is a square in, the more likely it is a square in the integers.

Use the PoC: QSv4_001.py -keysize 40

It is very slowly, but I want to keep working in this direction. I will continue that work here while finishing v3.
