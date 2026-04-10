NUMERICAL INTEGRATION CALCULATOR
Version: 1.0 (Midterm)

---

## OVERVIEW

This application is a desktop-based Numerical Integration Calculator built using Python and Tkinter. It computes definite integrals using multiple numerical methods and provides step-by-step solution trails, graphical visualization, and accuracy verification.

---

## FEATURES

* Trapezoidal Rule
* Simpson’s 1/3 Rule
* Simpson’s 3/8 Rule
* Midpoint Rule
* Step-by-step Solution Trail
* Graph Visualization of Functions
* Accuracy Verification:
  • Richardson Extrapolation
  • SciPy Adaptive Quadrature
  • SymPy Exact Solution
* Input Validation and Error Handling

---

## HOW TO RUN

1. Make sure Python 3.x is installed.
2. Install required libraries:
   pip install numpy scipy sympy matplotlib
3. Run the program:
   python main.py

---

## HOW TO USE

1. Enter a function (e.g., sin(x), x**2, exp(x))
2. Input lower bound (a) and upper bound (b)
3. Enter number of intervals (n)
4. Click "Calculate"
5. View results for each method
6. Click:

   * "View Solution" for step-by-step process
   * "Plot" for graphical visualization

---

## SUPPORTED FUNCTIONS

* sin, cos, tan
* exp, log, sqrt
* abs
* Constants: pi, e
* Variable: x

---

## EDGE CASE HANDLING

* Invalid function input → error message
* a ≥ b → prevented
* n < 2 → prevented
* Simpson constraints automatically adjusted
* Division by zero handled
* Domain errors handled (e.g., log(-1))

---

## LIMITATIONS

* Uses eval() for function parsing (not secure for production)
* Accuracy depends on number of intervals (n)
* Discontinuous functions may still produce approximation issues

---

## DEVELOPERS

* Gabriel Estrella
* Daniel Christian Mendoza
* Paul Rosal

---

## DESCRIPTION

This system is designed not only to compute numerical integrals but also to help users understand the process through solution trails and accuracy verification, making it both a computational and educational tool.

---

## END OF FILE
