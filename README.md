# Demo of small vibrations of a spherical membrane

This repo is for educational purposes.

### Introduction
The vibrations of a spherical membrane with radius $R$ are described by the following equation

$$v_{tt} = \bigtriangleup v\$$

in spherical coordinates. Let's consider the solutions that are in the form $v(r,\phi, \theta,t) \approx R + u(\phi, \theta, t)$ for $r \rightarrow R$, then

$$u_{tt} = c^2 \left(u_{\phi \phi} + \cot \phi u_\phi + \frac{1}{\sin^2 \phi} u_{\theta \theta} \right)$$

The general solution is
$$u(\theta, \phi, t) = \sum_n \sum_m a_{mm} u_{mn}(\theta, \phi, t)$$
where

$$u_{mn}(\theta, \phi, t) = P_n^m(\cos \phi) \cos(m\theta) \cos(\omega_{mn} t)$$

with $m \leq n$ and $\omega_{mn} = c \sqrt{n(n+1)}$ and $P_n^m$ are the Associated Legendre Polynomials.

The code calculates $(R + u_{mn})$ for $m=0,\dots,n_{max}$ and for $t\in[0,t_{max}]$, where $t_{max} = 2\pi/c$.

### Associated Legendre Polynomials
We are using normalized Associated Legendre Polynomials. The calculation is done via a recursion relation following https://arxiv.org/pdf/1410.1748

Thus, in this code we are actually calculating

$$\overline{P_n^m}(x) = \sqrt{\frac{(2n+1)(n-m)!}{2\pi(n+m)!}} P_n^m (x)$$

Hence

$$\int_{-1}^{1} \overline{P_n^m}(x) \overline{P_n^m}(x) dx = \frac{1}{\pi}$$

### How to run

To calculate $u_{mn}$ you have to
```Bash
cargo run NMAX
```
The solutions are written in text files at 
```Bash
data/{n}_{m}.txt
```
After that, you can plot these solutions using a Python script
```Bash
python scripts/view_data.py NMAX NT
```
This script will create a panel plotting each $u_{mn}$ for multiples times. Use $NT=0$ to plot a single snapshot. Each frame has also a small rotation in the viewer angle. The figures are stored at
```Bash
fig/fig_{t}.txt
```
Finally, you can generate a gif or mp4 with
```Bash
bash scripts/make_gif.sh
```

This is a sample of the resulting gif using NMAX = 6 and NT = 120
![alt text](./fig/out.gif)
References: 
https://faculty.fiu.edu/~meziani/Note13.pdf
https://arxiv.org/pdf/1410.1748