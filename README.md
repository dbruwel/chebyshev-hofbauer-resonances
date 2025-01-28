# DENR3901: Denison Research

This project was undertaken during the summer of 2024/2025 under the supervision of Dr. Caroline Wormell at the University of Sydney.

## Project Overview
A one-dimensional discrete-time dynamic system is defined by the equation: $$x_{t+1}=f(x_t)$$ s some function. If the points $x_t$ follow a distribution $x_t\sim\psi$, then the distribution of points at time $t+1$ is given by $\mathcal{L}\psi$, where $L$ is the "transfer operator."

The eigenfunctions of this operator are important because they represent distributions that evolve only by being scaled by an eigenvalue $\lambda$. If $|\lambda|<1$, these modes decay, while if $|\lambda|=1$, the modes are stable distributions that do not continue to evolve.

This project focuses on computational techniques for estimating the eigenvalues (or resonances) of the transfer operator, specifically for non-Markovian systems. The method uses Chebyshev approximations and Hofbauer extensions.