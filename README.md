# factorization_v4

For the current fastest version, look at v3 (which is also still a work in progress as I implement 2d sieving and experiment with lifting)... that one will easily factor above 200 bits. And once I port my findings to C and do all the remaining optimizations I hope to blaze past the presumed 110 digit hard ceiling for quadratic sieve - based algorithms.

This is however the cutting edge of my number theoretical research and the direction I should continue in, I've added a short introduction to this approach in my paper at chapter 7, so please refer to the paper to understand what is going on and how one possible approach to solve this may be a lattice basis reduction algorithm. Once I finish wrapping up v3 I will focus on this.

Use the PoC: QSv4_001.py -keysize 40

It is very slowly, but I want to keep working in this direction. I will continue that work here while finishing v3.

ps: The uploaded PoC is very bad and sloppy. It's a frankenstein PoC bc of just modifying previous iterations. I know that.. it's about the math though not the code.

Side note: I guess at this rate I will eventually have to rebrand my paper to: Factorization - A three year journey. lol. I don't know what it takes to get people to take me serious... I really have no choice but to keep going until the bitter end... it's the only way at this point.. especially since I can't find employment anymore or 0day buyers. literally no way out of my current life except succeeding at this.



