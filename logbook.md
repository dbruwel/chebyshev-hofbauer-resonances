# Logbook

### Initial details
**Goal:** workable code that takes in a (possibly nonlinear) tent map $f$, a tree depth, and some sort of grid size specification, and outputs a matrix approximation of the transfer operator, and perhaps some information about the intervals one is approximating it on.

**Nonlinear Tent Map:** have a function on $[0,1]$ that is uniformly expanding $(\inf |f'|>1)$ on a finite number of "pieces" that make up $[0,1]$, but not necessarily differentiable or continuous across the pieces. But increasing on $[0, 1/2]$ and decreasing on $[1/2,1]$ with $f(0)=f(1)=0$
 - **Uniformly expanding:** a function $f$ is uniformly expanding on some interval $I$ if $\inf|f'|>1$.
 - **Across the pieces:** this means at the boundaries, so each piece is well behaved (in $C^1$ and uniformly expanding), but each boundary can be discontinuous.
 - **Increasing / decreasing:** where the derivative is defined (at all points except a finite number of boundaries between pieces), the derivative is positive / negative.

The standard tent map is only considered a "nonlinear tent map" by this definition for the scale parameter greater than 1. We care about uniformly expanding because it is one of the key ways to get chaos (though there are situations where its not needed, e.g. the logistic map).

### Refined details
I think f is going to need to take in the "pieces" that the function is defined on in order to do the expansion well.

**Idea:**
 - **Chebyshev polynomials:** for any function $f$, we can write it as an expansion of Chebyshev polynomials, that is $f=a_0T_0(x)+a_1T_1(x)+...$. For some operator $L$, we can consider $Lf=a_0LT_0(x)+a_1LT_1(x)+...$. This means we only need to understand how $L$ acts on $T_n(x)$. We can further expand $LT_j(x)$ in terms of Chebyshev polynomials. We can write $\hat{L}_{ij}$ for the $i$th coefficient of $LT_j(x). Hence we can consider $Lf\approx \hat{L}\hat{f}$. We can further truncate after some number of Chebyshev terms. This gives us a way of approximating an operator by simply considering how it acts on the Chebyshev terms. Now suppose that $Lf=\lambda f$, i.e. $\lambda$ is an eigenvalue. We'd also have that $\lambda$ is "close" to an eigenvalue of $\hat{L}$. So the spectrum of $\hat{L}$ is an approximation for the spectrum of $L$.

