# Root Finding Project

## Project Description
This project implements a Python class called `RootFindingProblem` for solving equations of the form:

f(x) = 0

using classical numerical methods.

The project was developed without using any library root solvers, following the project instructions. Only standard Python and `cmath` are used.

---

## Implemented Methods

The class includes the following numerical methods:

1. Bisection Method
2. Fixed-Point Iteration
3. Newton's Method
4. Secant Method
5. False Position Method (Regula Falsi)
6. Horner's Method
7. Muller's Method

In addition, the internal structure also includes `_steffensen()` to match the required class structure.

---

## Algorithms and How They Work

### 1. Bisection Method
This method starts with an interval `[a, b]` where `f(a)` and `f(b)` have opposite signs.  
It repeatedly divides the interval into two halves and keeps the half containing the root.

### 2. Fixed-Point Iteration
This method rewrites the equation into the form:

x = g(x)

Starting from an initial guess `x0`, it repeatedly computes:

x_(n+1) = g(x_n)

until convergence.

### 3. Newton's Method
Newton's method uses the formula:

x_(n+1) = x_n - f(x_n) / f'(x_n)

It usually converges very quickly when the initial guess is good.

### 4. Secant Method
The secant method is similar to Newton's method but does not require the derivative.  
It uses two starting values and approximates the derivative using a secant line.

### 5. False Position Method
This method also starts with an interval `[a, b]` where the function changes sign.  
It finds the next approximation using the x-intercept of the secant line through the endpoints.

### 6. Horner's Method
Horner's method evaluates a polynomial efficiently using nested multiplication.

For example:

P(x) = x^3 - 6x^2 + 11x - 6

can be evaluated faster using Horner's scheme.

### 7. Muller's Method
Muller's method uses three starting points and a quadratic interpolation.  
It can find both real and complex roots.

---

## File Structure

```text
root-finding-project/
│── root_finding.py
│── examples.py
│── README.md
```
---

## How to Run the Examples
Open a terminal in the project folder and run:

python examples.py

This will execute several examples demonstrating the implemented methods.

---

## Example Using `solve()`

```python
from root_finding import RootFindingProblem

f = lambda x: x**3 - x - 2
df = lambda x: 3*x**2 - 1

p = RootFindingProblem(f=f, df=df)

root = p.solve("newton", x0=1.5)

print("Root:", root)

The solve() method selects the appropriate algorithm based on the method name and returns the computed root.