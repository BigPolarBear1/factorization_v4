# factorization_v4

For the current fastest version and my paper, look at v3 (which is also still a work in progress as I rewrite that project to C)... that one will easily factor above 200 bits. And once I port my findings to C and do all the remaining optimizations I hope to blaze past the presumed 110 digit hard ceiling for quadratic sieve - based algorithms.

This is however the cutting edge of my number theoretical research and the direction I should continue in.

The current PoC creates a hashmap that shows for every linear coefficient mod p at which quadratic coefficient mod p it is found.
I.e if we have as valid solution 6x<sup>2</sup>+4x = 0 mod p then we will find 6 as a valid solution as quadratic coefficient for linear coefficient 4 in EVERY MOD P (and its powers), without exception.
This "idea" is similar to using legendre symbols to determine squaredness of an integers by checking if its a square mod multiple primes and the more primes moduli it is a square in, the more likely it is a square in the integers.

There is two approaches to creating the hashmap:

1. Key by linear coefficient
2. Key by quadratic coefficient

For number 1, you would then look for a small quadratic coefficient that appears for every mod p.
For number 2, you would then look for a small linear coefficient that appears for every mod p.

The main problem that remains is, can we construct an efficient algorithm to do this? I will keep grindng this direction

Use the PoC: QSv4_001.py -keysize 40

It is very slowly, but I want to keep working in this direction. I will continue that work here while finishing v3.

ps: The uploaded PoC is very bad and sloppy. It's a frankenstein PoC bc of just modifying previous iterations. I know that.. it's about the math though not the code.

pss: Yes I know it is very slow, but it doesn't matter, that's just how this type of research works. Small steps while figuring things out. I think the representation in the hashmap is good, I just need to figure out how to pull solutions from it efficiently.. that's the main bulk that requires further research. Anyway, first for the coming days I am going to focus on v3 and porting to C.

UPDATE: Actuaaaally, I was just let my mind wander, thinking back about lattice basis reduction. If you just limit the scope to just one linear coefficient pairing, and to factor N, the problem then becomes finding a quadratic coefficient such that we have a square in every mod p. But since the quadratic coefficient is just how many times N we add or subtract after squaring the linear coefficient... determining that amount it needs to be squared may be a good candidate for lattice basis reduction. Let me wrap up v3 and start trying this.

