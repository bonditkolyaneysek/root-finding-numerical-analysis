from root_finding import RootFindingProblem


def main():
    print("===== ROOT FINDING EXAMPLES =====\n")

    # -------------------------------------------------
    # 1. Bisection Method
    # f(x) = x^3 - x - 2, root near x = 1.521...
    # -------------------------------------------------
    f1 = lambda x: x**3 - x - 2
    p1 = RootFindingProblem(f=f1)
    root_bisection = p1.solve("bisection", a=1, b=2)
    print("1) Bisection Method")
    print("Root =", root_bisection)
    print()

    # -------------------------------------------------
    # 2. Fixed-Point Iteration
    # x = cos(x)
    # so g(x) = cos(x)
    # -------------------------------------------------
    import math
    g = lambda x: math.cos(x)
    p2 = RootFindingProblem(g=g)
    root_fixed = p2.solve("fixed-point", x0=0.5)
    print("2) Fixed-Point Iteration")
    print("Root =", root_fixed)
    print()

    # -------------------------------------------------
    # 3. Newton's Method
    # f(x) = x^3 - x - 2
    # df(x) = 3x^2 - 1
    # -------------------------------------------------
    df1 = lambda x: 3 * x**2 - 1
    p3 = RootFindingProblem(f=f1, df=df1)
    root_newton = p3.solve("newton", x0=1.5)
    print("3) Newton's Method")
    print("Root =", root_newton)
    print()

    # -------------------------------------------------
    # 4. Secant Method
    # same function
    # -------------------------------------------------
    p4 = RootFindingProblem(f=f1)
    root_secant = p4.solve("secant", x0=1, x1=2)
    print("4) Secant Method")
    print("Root =", root_secant)
    print()

    # -------------------------------------------------
    # 5. False Position Method
    # same function
    # -------------------------------------------------
    p5 = RootFindingProblem(f=f1)
    root_false_position = p5.solve("false_position", a=1, b=2)
    print("5) False Position Method")
    print("Root =", root_false_position)
    print()

    # -------------------------------------------------
    # 6. Horner's Method
    # Evaluate polynomial x^3 - 6x^2 + 11x - 6 at x = 2
    # Expected result = 0
    # -------------------------------------------------
    p6 = RootFindingProblem()
    coeffs = [1, -6, 11, -6]
    horner_value = p6.solve("horner", coeffs=coeffs, x=2)
    print("6) Horner's Method")
    print("P(2) =", horner_value)
    print()

    # -------------------------------------------------
    # 7. Muller's Method
    # f(x) = x^2 + 1 has roots i and -i
    # -------------------------------------------------
    f_complex = lambda x: x**2 + 1
    p7 = RootFindingProblem(f=f_complex)
    root_muller = p7.solve("muller", x0=0, x1=1, x2=1 + 1j)
    print("7) Muller's Method")
    print("Complex root =", root_muller)
    print()

    print("All examples ran successfully.")


if __name__ == "__main__":
    main()