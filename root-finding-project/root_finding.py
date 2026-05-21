import cmath


class RootFindingProblem:
    def __init__(self, f=None, df=None, g=None):
        """
        f(x): function whose root is sought
        df(x): derivative of f(x), used in Newton's method
        g(x): fixed-point function, used in fixed-point iteration
        """
        self.f = f
        self.df = df
        self.g = g

    def solve(self, method, **kwargs):
        """
        Public solver interface.
        Example:
            p.solve("bisection", a=0, b=1)
            p.solve("newton", x0=0.5)
            p.solve("secant", x0=0, x1=1)
        """
        method = method.lower().strip()

        if method == "bisection":
            return self._bisection(**kwargs)
        elif method in ("fixed-point", "fixed_point", "fixedpoint"):
            return self._fixed_point(**kwargs)
        elif method == "newton":
            return self._newton(**kwargs)
        elif method == "secant":
            return self._secant(**kwargs)
        elif method in ("false-position", "false_position", "regula_falsi", "regular_falsi"):
            return self._false_position(**kwargs)
        elif method == "horner":
            return self._horner(**kwargs)
        elif method == "muller":
            return self._muller(**kwargs)
        else:
            raise ValueError(f"Unknown method: {method}")

    def _bisection(self, a, b, tol=1e-10, max_iter=100):
        if self.f is None:
            raise ValueError("Function f is required for bisection method.")

        fa = self.f(a)
        fb = self.f(b)

        if fa * fb > 0:
            raise ValueError("Invalid interval: f(a) and f(b) must have opposite signs.")

        for _ in range(max_iter):
            c = (a + b) / 2
            fc = self.f(c)

            if abs(fc) < tol or abs(b - a) / 2 < tol:
                return c

            if fa * fc < 0:
                b = c
                fb = fc
            else:
                a = c
                fa = fc

        raise ValueError("No convergence after max_iter in bisection method.")

    def _fixed_point(self, x0, tol=1e-10, max_iter=100):
        if self.g is None:
            raise ValueError("Missing fixed-point function g for fixed-point iteration.")

        x = x0
        for _ in range(max_iter):
            x_next = self.g(x)

            if abs(x_next - x) < tol:
                return x_next

            x = x_next

        raise ValueError("No convergence after max_iter in fixed-point iteration.")

    def _newton(self, x0, tol=1e-10, max_iter=100):
        if self.f is None:
            raise ValueError("Function f is required for Newton's method.")
        if self.df is None:
            raise ValueError("Missing derivative df for Newton's method.")

        x = x0
        for _ in range(max_iter):
            fx = self.f(x)
            dfx = self.df(x)

            if dfx == 0:
                raise ZeroDivisionError("Division by zero during Newton iteration.")

            x_next = x - fx / dfx

            if abs(x_next - x) < tol or abs(self.f(x_next)) < tol:
                return x_next

            x = x_next

        raise ValueError("No convergence after max_iter in Newton's method.")

    def _secant(self, x0, x1, tol=1e-10, max_iter=100):
        if self.f is None:
            raise ValueError("Function f is required for secant method.")

        for _ in range(max_iter):
            f0 = self.f(x0)
            f1 = self.f(x1)

            if f1 - f0 == 0:
                raise ZeroDivisionError("Division by zero during secant iteration.")

            x2 = x1 - f1 * (x1 - x0) / (f1 - f0)

            if abs(x2 - x1) < tol or abs(self.f(x2)) < tol:
                return x2

            x0, x1 = x1, x2

        raise ValueError("No convergence after max_iter in secant method.")

    def _false_position(self, a, b, tol=1e-10, max_iter=100):
        if self.f is None:
            raise ValueError("Function f is required for false position method.")

        fa = self.f(a)
        fb = self.f(b)

        if fa * fb > 0:
            raise ValueError("Invalid interval: f(a) and f(b) must have opposite signs.")

        for _ in range(max_iter):
            denominator = fb - fa
            if denominator == 0:
                raise ZeroDivisionError("Division by zero during false position iteration.")

            c = (a * fb - b * fa) / (fb - fa)
            fc = self.f(c)

            if abs(fc) < tol or abs(b - a) < tol:
                return c

            if fa * fc < 0:
                b = c
                fb = fc
            else:
                a = c
                fa = fc

        raise ValueError("No convergence after max_iter in false position method.")

    def _steffensen(self, x0, tol=1e-10, max_iter=100):
        """
        Included because it appears in the required class structure in the PDF.
        Not part of the listed required public methods, but implemented to match
        the requested internal structure.
        """
        if self.g is None:
            raise ValueError("Missing function g for Steffensen's method.")

        x = x0
        for _ in range(max_iter):
            gx = self.g(x)
            ggx = self.g(gx)
            denominator = ggx - 2 * gx + x

            if denominator == 0:
                raise ZeroDivisionError("Division by zero during Steffensen iteration.")

            x_next = x - ((gx - x) ** 2) / denominator

            if abs(x_next - x) < tol:
                return x_next

            x = x_next

        raise ValueError("No convergence after max_iter in Steffensen's method.")

    def _horner(self, coeffs, x):
        """
        Evaluate polynomial using Horner's method.
        coeffs example for x^3 - 6x^2 + 11x - 6:
            [1, -6, 11, -6]
        """
        if not coeffs:
            raise ValueError("Coefficient list cannot be empty.")

        result = coeffs[0]
        for coeff in coeffs[1:]:
            result = result * x + coeff
        return result

    def _muller(self, x0, x1, x2, tol=1e-10, max_iter=100):
        if self.f is None:
            raise ValueError("Function f is required for Muller's method.")

        x0 = complex(x0)
        x1 = complex(x1)
        x2 = complex(x2)

        for _ in range(max_iter):
            h1 = x1 - x0
            h2 = x2 - x1

            if h1 == 0 or h2 == 0:
                raise ZeroDivisionError("Division by zero during Muller iteration.")

            d1 = (self.f(x1) - self.f(x0)) / h1
            d2 = (self.f(x2) - self.f(x1)) / h2

            denominator_h = h2 + h1
            if denominator_h == 0:
                raise ZeroDivisionError("Division by zero during Muller iteration.")

            d = (d2 - d1) / denominator_h
            b = d2 + h2 * d
            c = self.f(x2)

            rad = cmath.sqrt(b * b - 4 * c * d)

            if abs(b + rad) > abs(b - rad):
                den = b + rad
            else:
                den = b - rad

            if den == 0:
                raise ZeroDivisionError("Division by zero during Muller iteration.")

            h = -2 * c / den
            x3 = x2 + h

            if abs(h) < tol or abs(self.f(x3)) < tol:
                return x3

            x0, x1, x2 = x1, x2, x3

        raise ValueError("No convergence after max_iter in Muller's method.")