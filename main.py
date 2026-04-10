
import tkinter as tk
from tkinter import messagebox
import math
import traceback
import numpy as np
from scipy import integrate as sci_integrate
import sympy as sp
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.patches as patches
import matplotlib.colors as mcolors

# ── Colors ───────────────────────────────────────────────────────
BG     = "#EEEEF8"
PANEL  = "#FFFFFF"
PURPLE = "#5B5BD6"
GREEN  = "#2B7A0B"
ORANGE = "#C75C00"
VIOLET = "#7C66DC"
TEXT   = "#1C1C2E"
DIM    = "#7070A0"
RED    = "#CC2222"
BORDER = "#DDDDF0"

METHODS = [
    ("Trapezoidal Rule",   "#5B5BD6", "Uses trapezoids for approximation"),
    ("Simpson's 1/3 Rule", "#7C66DC", "Uses parabolic segments"),
    ("Midpoint Rule",      "#2B7A0B", "Uses rectangle midpoints"),
    ("Simpson's 3/8 Rule", "#C75C00", "Uses cubic polynomials"),
]

# ── Math ─────────────────────────────────────────────────────────
NS = {k: v for k, v in vars(math).items() if not k.startswith("_")}
NS.update({"__builtins__": {}, "e": math.e, "pi": math.pi})

def make_f(expr):
    def f(x):
        env = dict(NS); env["x"] = x
        return float(eval(expr, env))
    return f

def clean(raw):
    return raw.strip().replace("^", "**")

def trap(f, a, b, n):
    h = (b-a)/n
    xs = [a + i*h for i in range(n+1)]
    ys = [f(x) for x in xs]
    return h * (ys[0]/2 + sum(ys[1:-1]) + ys[-1]/2)

def s13(f, a, b, n):
    if n % 2: n += 1
    h = (b-a)/n
    xs = [a + i*h for i in range(n+1)]
    ys = [f(x) for x in xs]
    coeffs = [1] + [4 if i%2 else 2 for i in range(1, n)] + [1]
    return h/3 * sum(c*y for c,y in zip(coeffs, ys))

def s38(f, a, b, n):
    if n % 3: n = (n//3+1)*3
    h = (b-a)/n
    xs = [a + i*h for i in range(n+1)]
    ys = [f(x) for x in xs]
    coeffs = [1] + [3 if i%3 else 2 for i in range(1, n)] + [1]
    return 3*h/8 * sum(c*y for c,y in zip(coeffs, ys))

def mid(f, a, b, n):
    h = (b-a)/n
    return h * sum(f(a + (i+0.5)*h) for i in range(n))

FUNCS = {
    "Trapezoidal Rule":   trap,
    "Simpson's 1/3 Rule": s13,
    "Midpoint Rule":      mid,
    "Simpson's 3/8 Rule": s38,
}

# ── Solution Trail builders ───────────────────────────────────────
def make_trail_trap(raw, expr, a, b, n, result):
    f = make_f(expr)
    h = (b-a)/n
    xs = [a + i*h for i in range(n+1)]
    ys = [f(x) for x in xs]
    int_sum = sum(ys[1:-1])
    weighted = ys[0] + 2*int_sum + ys[-1]

    L = []
    L.append(("h", "TRAPEZOIDAL RULE - SOLUTION TRAIL"))
    L.append(("", ""))
    L.append(("s", "GIVEN"))
    L.append(("b", f"  f(x) = {raw}"))
    L.append(("b", f"  a = {a},  b = {b},  n = {n}"))
    L.append(("", ""))
    L.append(("s", "STEP 1 - Interval width h"))
    L.append(("b", f"  h = (b - a) / n  =  ({b} - {a}) / {n}"))
    L.append(("r", f"  h = {h:.8f}"))
    L.append(("", ""))
    L.append(("s", "STEP 2 - Formula"))
    L.append(("b", "  Integral ~ (h/2) x [f(x0) + 2f(x1) + ... + 2f(xn-1) + f(xn)]"))
    L.append(("b", "  Endpoint coefficients = 1,  Interior = 2"))
    L.append(("", ""))
    L.append(("s", "STEP 3 - Function evaluations"))
    L.append(("b", f"  f(x0) = f({a}) = {ys[0]:.8f}   [coeff=1]"))
    for i in range(1, min(n, 8)+1):
        L.append(("b", f"  f(x{i}) = f({xs[i]:.5f}) = {ys[i]:.8f}   [coeff=2]"))
    if n > 8:
        L.append(("d", f"  ... {n-8} more interior points ..."))
    L.append(("b", f"  f(x{n}) = f({b}) = {ys[-1]:.8f}   [coeff=1]"))
    L.append(("", ""))
    L.append(("s", "STEP 4 - Weighted sum"))
    L.append(("b", f"  Sum = {ys[0]:.8f} + 2 x {int_sum:.8f} + {ys[-1]:.8f}"))
    L.append(("r", f"  Sum = {weighted:.8f}"))
    L.append(("", ""))
    L.append(("s", "STEP 5 - Multiply by h/2"))
    L.append(("b", f"  Result = (h/2) x Sum  =  ({h:.6f}/2) x {weighted:.6f}"))
    L.append(("r", f"  Result = {result:.10f}"))
    L.append(("", ""))
    _append_verify(L, trap, expr, raw, a, b, n, result)
    return L

def make_trail_s13(raw, expr, a, b, n, result):
    if n % 2: n += 1
    f = make_f(expr)
    h = (b-a)/n
    xs = [a + i*h for i in range(n+1)]
    ys = [f(x) for x in xs]
    coeffs = [1] + [4 if i%2 else 2 for i in range(1, n)] + [1]
    weighted = sum(c*y for c,y in zip(coeffs, ys))

    L = []
    L.append(("h", "SIMPSON'S 1/3 RULE - SOLUTION TRAIL"))
    L.append(("", ""))
    L.append(("s", "GIVEN"))
    L.append(("b", f"  f(x) = {raw}"))
    L.append(("b", f"  a = {a},  b = {b},  n = {n}  (must be even)"))
    L.append(("", ""))
    L.append(("s", "STEP 1 - Verify n is even"))
    L.append(("b", f"  n = {n}  -> {'OK (even)' if n%2==0 else 'adjusted'}"))
    L.append(("", ""))
    L.append(("s", "STEP 2 - Interval width h"))
    L.append(("r", f"  h = ({b} - {a}) / {n} = {h:.8f}"))
    L.append(("", ""))
    L.append(("s", "STEP 3 - Formula"))
    L.append(("b", "  Integral ~ (h/3) x [f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + ... + f(xn)]"))
    L.append(("b", "  Pattern: 1, 4, 2, 4, 2, ..., 4, 1"))
    L.append(("", ""))
    L.append(("s", "STEP 4 - Evaluations with coefficients"))
    for i in range(min(n+1, 9)):
        L.append(("b", f"  f(x{i}) = f({xs[i]:.5f}) = {ys[i]:.8f}   [coeff={coeffs[i]}]"))
    if n+1 > 9:
        L.append(("d", f"  ... {n+1-9} more points ..."))
    L.append(("", ""))
    L.append(("s", "STEP 5 - Weighted sum"))
    L.append(("r", f"  Weighted sum = {weighted:.8f}"))
    L.append(("", ""))
    L.append(("s", "STEP 6 - Multiply by h/3"))
    L.append(("b", f"  Result = (h/3) x {weighted:.6f}  =  ({h:.6f}/3) x {weighted:.6f}"))
    L.append(("r", f"  Result = {result:.10f}"))
    L.append(("", ""))
    _append_verify(L, s13, expr, raw, a, b, n, result)
    return L

def make_trail_mid(raw, expr, a, b, n, result):
    f = make_f(expr)
    h = (b-a)/n
    mids = [a + (i+0.5)*h for i in range(n)]
    ys   = [f(m) for m in mids]

    L = []
    L.append(("h", "MIDPOINT RULE - SOLUTION TRAIL"))
    L.append(("", ""))
    L.append(("s", "GIVEN"))
    L.append(("b", f"  f(x) = {raw}"))
    L.append(("b", f"  a = {a},  b = {b},  n = {n}"))
    L.append(("", ""))
    L.append(("s", "STEP 1 - Interval width h"))
    L.append(("r", f"  h = ({b} - {a}) / {n} = {h:.8f}"))
    L.append(("", ""))
    L.append(("s", "STEP 2 - Formula"))
    L.append(("b", "  Integral ~ h x SUM of f(midpoints)"))
    L.append(("b", "  Midpoint of interval i = a + (i + 0.5) x h"))
    L.append(("", ""))
    L.append(("s", "STEP 3 - Midpoint evaluations"))
    for i in range(min(n, 8)):
        L.append(("b", f"  Interval {i}: mid={mids[i]:.6f}  f={ys[i]:.8f}"))
    if n > 8:
        L.append(("d", f"  ... {n-8} more intervals ..."))
    L.append(("", ""))
    L.append(("s", "STEP 4 - Sum"))
    L.append(("r", f"  SUM = {sum(ys):.8f}"))
    L.append(("", ""))
    L.append(("s", "STEP 5 - Multiply by h"))
    L.append(("b", f"  Result = {h:.8f} x {sum(ys):.8f}"))
    L.append(("r", f"  Result = {result:.10f}"))
    L.append(("", ""))
    _append_verify(L, mid, expr, raw, a, b, n, result)
    return L

def make_trail_s38(raw, expr, a, b, n, result):
    if n % 3: n = (n//3+1)*3
    f = make_f(expr)
    h = (b-a)/n
    xs = [a + i*h for i in range(n+1)]
    ys = [f(x) for x in xs]
    coeffs = [1] + [3 if i%3 else 2 for i in range(1, n)] + [1]
    weighted = sum(c*y for c,y in zip(coeffs, ys))

    L = []
    L.append(("h", "SIMPSON'S 3/8 RULE - SOLUTION TRAIL"))
    L.append(("", ""))
    L.append(("s", "GIVEN"))
    L.append(("b", f"  f(x) = {raw}"))
    L.append(("b", f"  a = {a},  b = {b},  n = {n}  (multiple of 3)"))
    L.append(("", ""))
    L.append(("s", "STEP 1 - Interval width h"))
    L.append(("r", f"  h = ({b} - {a}) / {n} = {h:.8f}"))
    L.append(("", ""))
    L.append(("s", "STEP 2 - Formula"))
    L.append(("b", "  Integral ~ (3h/8) x [f(x0)+3f(x1)+3f(x2)+2f(x3)+...+f(xn)]"))
    L.append(("b", "  Pattern: 1, 3, 3, 2, 3, 3, 2, ..., 1"))
    L.append(("", ""))
    L.append(("s", "STEP 3 - Evaluations with coefficients"))
    for i in range(min(n+1, 9)):
        L.append(("b", f"  f(x{i}) = f({xs[i]:.5f}) = {ys[i]:.8f}   [coeff={coeffs[i]}]"))
    if n+1 > 9:
        L.append(("d", f"  ... {n+1-9} more points ..."))
    L.append(("", ""))
    L.append(("s", "STEP 4 - Result"))
    L.append(("b", f"  Weighted sum = {weighted:.8f}"))
    L.append(("b", f"  Result = (3 x {h:.6f} / 8) x {weighted:.6f}"))
    L.append(("r", f"  Result = {result:.10f}"))
    L.append(("", ""))
    _append_verify(L, s38, expr, raw, a, b, n, result)
    return L

def _append_verify(L, fn, expr, raw, a, b, n, result):
    f = make_f(expr)
    r2 = fn(f, a, b, n*2)
    ref, _ = sci_integrate.quad(f, a, b)
    abs_err = abs(result - ref)
    rich_err = abs(result - r2)

    L.append(("sep", "-" * 52))
    L.append(("v",   "VERIFICATION & ACCURACY AUDIT"))
    L.append(("", ""))
    L.append(("b",   "  [1] Richardson Extrapolation:"))
    L.append(("b",   f"      n={n}: {result:.10f}"))
    L.append(("b",   f"      n={n*2}: {r2:.10f}"))
    L.append(("b",   f"      Estimated error = {rich_err:.4e}"))
    L.append(("", ""))
    L.append(("b",   "  [2] SciPy adaptive quadrature:"))
    L.append(("b",   f"      Reference = {ref:.10f}"))
    L.append(("b",   f"      Absolute error = {abs_err:.4e}"))
    rel = abs_err/abs(ref)*100 if ref != 0 else 0
    L.append(("b",   f"      Relative error = {rel:.6f}%"))
    L.append(("", ""))

    try:
        x = sp.Symbol("x")
        ex = sp.sympify(raw.replace("^","**"))
        ai = sp.integrate(ex, x)
        exact = float(ai.subs(x, b) - ai.subs(x, a))
        L.append(("b", "  [3] SymPy exact symbolic:"))
        L.append(("b", f"      Antiderivative = {ai}"))
        L.append(("b", f"      Exact value    = {exact:.10f}"))
        L.append(("b", f"      Error vs exact = {abs(result-exact):.4e}"))
        L.append(("", ""))
    except Exception:
        pass

    ok = abs_err < 1e-4
    L.append(("ok" if ok else "warn",
               "  RESULT: HIGH ACCURACY" if ok else "  RESULT: CONSIDER INCREASING n"))
    L.append(("", ""))
    L.append(("sep", "-" * 52))
    L.append(("h",   "SUMMARY"))
    L.append(("r",   f"  Integral from {a} to {b} of ({raw}) dx"))
    L.append(("r",   f"  = {result:.10f}"))
    L.append(("b",   f"  Error ~ {abs_err:.4e}"))

TRAIL_BUILDERS = {
    "Trapezoidal Rule":   make_trail_trap,
    "Simpson's 1/3 Rule": make_trail_s13,
    "Midpoint Rule":      make_trail_mid,
    "Simpson's 3/8 Rule": make_trail_s38,
}

# ── Solution Trail Window ─────────────────────────────────────────
class TrailWindow(tk.Toplevel):
    def __init__(self, parent, method, color, lines):
        super().__init__(parent)
        self.title(f"Solution Trail - {method}")
        self.geometry("640x600")
        self.resizable(True, True)
        self.configure(bg=color)
        self.grab_set()

        # Header
        hdr = tk.Frame(self, bg=color)
        hdr.pack(fill="x", ipady=8)
        tk.Label(hdr, text="Solution Trail", bg=color, fg="white",
                 font=("Segoe UI", 13, "bold")).pack(side="left", padx=16)
        tk.Label(hdr, text=method, bg=color, fg="#CCCCEE",
                 font=("Segoe UI", 9)).pack(side="left")
        tk.Button(hdr, text="Close", bg=color, fg="white",
                  relief="flat", font=("Segoe UI", 9),
                  cursor="hand2", command=self.destroy
                  ).pack(side="right", padx=12)

        # Text area
        frm = tk.Frame(self, bg="white")
        frm.pack(fill="both", expand=True)

        vsb = tk.Scrollbar(frm)
        vsb.pack(side="right", fill="y")
        hsb = tk.Scrollbar(frm, orient="horizontal")
        hsb.pack(side="bottom", fill="x")

        txt = tk.Text(frm, wrap="none",
                      yscrollcommand=vsb.set, xscrollcommand=hsb.set,
                      bg="#F9F9FF", fg=TEXT,
                      font=("Courier New", 10),
                      relief="flat", padx=18, pady=14,
                      borderwidth=0)
        txt.pack(fill="both", expand=True)
        vsb.config(command=txt.yview)
        hsb.config(command=txt.xview)
        txt.bind("<MouseWheel>",
                 lambda e: txt.yview_scroll(-1*int(e.delta/120), "units"))

        # Tags
        txt.tag_config("h",    foreground=color,       font=("Courier New", 10, "bold"))
        txt.tag_config("s",    foreground=color,       font=("Courier New", 10, "bold"))
        txt.tag_config("v",    foreground=GREEN,       font=("Courier New", 10, "bold"))
        txt.tag_config("r",    foreground="#111133",   font=("Courier New", 10, "bold"))
        txt.tag_config("b",    foreground="#222244",   font=("Courier New", 10))
        txt.tag_config("d",    foreground="#9999BB",   font=("Courier New", 10))
        txt.tag_config("sep",  foreground="#CCCCDD",   font=("Courier New", 10))
        txt.tag_config("ok",   foreground=GREEN,       font=("Courier New", 10, "bold"))
        txt.tag_config("warn", foreground=ORANGE,      font=("Courier New", 10, "bold"))
        txt.tag_config("",     foreground="#F9F9FF",   font=("Courier New", 4))

        for tag, text in lines:
            txt.insert("end", text + "\n", tag)
        txt.config(state="disabled")


# ── Chart Window ──────────────────────────────────────────────────
class ChartWindow(tk.Toplevel):
    def __init__(self, parent, raw, expr, a, b, n, method, color):
        super().__init__(parent)
        self.title(f"Plot - {method}")
        self.geometry("720x480")
        self.resizable(True, True)
        f   = make_f(expr)
        res = FUNCS[method](f, a, b, n)
        rgb = mcolors.to_rgb(color)
        h   = (b - a) / n

        fig = Figure(figsize=(7, 4.4), dpi=100, facecolor="#FAFAFA")
        ax  = fig.add_subplot(111, facecolor="#F5F5FF")
        xs  = np.linspace(a, b, 500)
        ys  = np.vectorize(f)(xs)

        if method == "Midpoint Rule":
            for i in range(min(n, 60)):
                xi, xi1 = a+i*h, a+(i+1)*h
                ht = f((xi+xi1)/2)
                ax.add_patch(patches.Rectangle(
                    (xi, min(0,ht)), h, abs(ht),
                    lw=0.5, edgecolor=color, facecolor=(*rgb, 0.2)))
        elif method == "Trapezoidal Rule":
            xs_t = np.linspace(a, b, min(n,60)+1)
            ys_t = np.vectorize(f)(xs_t)
            for i in range(min(n, 60)):
                ax.fill([xs_t[i],xs_t[i+1],xs_t[i+1],xs_t[i]],
                        [0,0,ys_t[i+1],ys_t[i]], color=(*rgb, 0.2), lw=0)
                ax.plot([xs_t[i],xs_t[i+1]], [ys_t[i],ys_t[i+1]], color=color, lw=0.5)
        else:
            ax.fill_between(xs, ys, alpha=0.18, color=color)

        ax.plot(xs, ys, color=color, lw=2, label=f"f(x) = {raw}")
        ax.axhline(0, color="#AAAACC", lw=0.8)
        ax.set_title(f"Integral from {a} to {b} of {raw}  ~  {res:.6f}", fontsize=10)
        ax.legend(fontsize=9); ax.spines[["top","right"]].set_visible(False)
        fig.tight_layout(pad=1.4)
        FigureCanvasTkAgg(fig, master=self).get_tk_widget().pack(fill="both", expand=True)
        FigureCanvasTkAgg(fig, master=self).draw()

class AboutWindow(tk.Toplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.title("About / Help")
        self.geometry("520x500")
        self.configure(bg="#FFFFFF")
        self.grab_set()

        # ── Scrollable Canvas ─────────────────────────
        canvas = tk.Canvas(self, bg="#FFFFFF", highlightthickness=0)
        scrollbar = tk.Scrollbar(self, orient="vertical", command=canvas.yview)
        scroll_frame = tk.Frame(canvas, bg="#FFFFFF")

        scroll_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # ── Content Container ─────────────────────────
        container = tk.Frame(scroll_frame, bg="#FFFFFF", padx=20, pady=20)
        container.pack(fill="both", expand=True)

        # Title
        tk.Label(container,
                 text="Numerical Integration Calculator",
                 font=("Segoe UI", 15, "bold"),
                 bg="#FFFFFF", fg="#1C1C2E").pack(anchor="w", pady=(0,10))

        # Version
        tk.Label(container,
                 text="Version: 1.0 (Midterm)",
                 font=("Segoe UI", 10),
                 bg="#FFFFFF", fg="#555577").pack(anchor="w", pady=(0,10))

        # Description
        tk.Label(container,
                 text="This application computes definite integrals using "
                      "numerical methods and provides step-by-step solutions "
                      "and accuracy verification.",
                 wraplength=460,
                 justify="left",
                 font=("Segoe UI", 10),
                 bg="#FFFFFF", fg="#333344").pack(anchor="w", pady=(0,15))

        # Features
        tk.Label(container, text="Features:",
                 font=("Segoe UI", 11, "bold"),
                 bg="#FFFFFF").pack(anchor="w")

        features = [
            "Trapezoidal Rule",
            "Simpson’s 1/3 Rule",
            "Simpson’s 3/8 Rule",
            "Midpoint Rule",
            "Solution Trail (step-by-step)",
            "Graph Visualization",
            "Accuracy Verification (SciPy, SymPy)"
        ]

        for f in features:
            tk.Label(container, text=f"• {f}",
                     font=("Segoe UI", 10),
                     bg="#FFFFFF", fg="#444466").pack(anchor="w")

        # Members
        tk.Label(container, text="\nDeveloped by:",
                 font=("Segoe UI", 11, "bold"),
                 bg="#FFFFFF").pack(anchor="w")

        tk.Label(container,
                 text="• Gabriel Estrella\n• Daniel Christian Mendoza\n• Paul Rosal",
                 justify="left",
                 font=("Segoe UI", 10),
                 bg="#FFFFFF").pack(anchor="w")

        # How to Use
        tk.Label(container, text="\nHow to Use:",
                 font=("Segoe UI", 11, "bold"),
                 bg="#FFFFFF").pack(anchor="w")

        tk.Label(container,
                 text="1. Enter a function (e.g., sin(x), x**2)\n"
                      "2. Set lower and upper bounds\n"
                      "3. Enter number of intervals\n"
                      "4. Click Calculate\n"
                      "5. View solution or plot",
                 justify="left",
                 wraplength=460,
                 font=("Segoe UI", 10),
                 bg="#FFFFFF").pack(anchor="w")

        # Close button
        tk.Button(container, text="Close",
                  bg="#5B5BD6", fg="white",
                  font=("Segoe UI", 10, "bold"),
                  relief="flat",
                  command=self.destroy).pack(pady=15)

        # ── Mouse Scroll Support ───────────────────────
        canvas.bind_all("<MouseWheel>",
            lambda e: canvas.yview_scroll(-1 * int(e.delta / 120), "units"))
        s
# ── Main App ──────────────────────────────────────────────────────
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Numerical Integration Calculator")
        self.geometry("1050x700")
        self.configure(bg=BG)
        self.resizable(True, True)
        self._results  = {}
        self._computed = {}
        self._build()

    def _build(self):
        # ── Title bar ─────────────────────────────────────
        top = tk.Frame(self, bg=BG)
        top.pack(fill="x", padx=20, pady=14)
        icon = tk.Frame(top, bg=PURPLE, width=48, height=48)
        icon.pack(side="left"); icon.pack_propagate(False)
        tk.Label(icon, text="∫", bg=PURPLE, fg="white",
                 font=("Segoe UI",24,"bold")).place(relx=.5,rely=.5,anchor="center")
        tf = tk.Frame(top, bg=BG)
        tf.pack(side="left", padx=10)
        tk.Label(tf, text="Numerical Integration Calculator",
                 bg=BG, fg=TEXT, font=("Segoe UI",17,"bold")).pack(anchor="w")
        tk.Label(tf, text="Calculate definite integrals with multiple methods",
                 bg=BG, fg=DIM, font=("Segoe UI",9)).pack(anchor="w")
        
        tk.Button(top,
          text="About / Help",
          bg="#EEEEF8",
          fg="#5B5BD6",
          font=("Segoe UI", 9, "bold"),
          relief="flat",
          cursor="hand2",
          command=lambda: AboutWindow(self)
          ).pack(side="right", padx=10)

        # ── Two columns ───────────────────────────────────
        body = tk.Frame(self, bg=BG)
        body.pack(fill="both", expand=True, padx=18, pady=14)

        # LEFT panel (fixed width, scrollable)
        left_frame = tk.Frame(body, bg=BG, width=300)
        left_frame.pack(side="left", fill="y", padx=14)
        left_frame.pack_propagate(False)

        canvas = tk.Canvas(left_frame, bg=BG, highlightthickness=0, width=284)
        vsb    = tk.Scrollbar(left_frame, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        left = tk.Frame(canvas, bg=BG)
        win  = canvas.create_window((0,0), window=left, anchor="nw")

        def _resize(e):
            canvas.configure(scrollregion=canvas.bbox("all"))
            canvas.itemconfig(win, width=e.width if hasattr(e,'width') else canvas.winfo_width())

        left.bind("<Configure>", _resize)
        canvas.bind("<Configure>", lambda e: canvas.itemconfig(win, width=e.width))
        canvas.bind_all("<MouseWheel>",
                        lambda e: canvas.yview_scroll(-1*int(e.delta/120), "units"))

        # RIGHT panel
        right = tk.Frame(body, bg=BG)
        right.pack(side="left", fill="both", expand=True)

        self._build_left(left)
        self._build_right(right)

    # ── Left panel ────────────────────────────────────────
    def _build_left(self, parent):
        # Setup card
        card = self._card(parent)
        card.pack(fill="x", pady=10, padx=2)
        self._hdr(card, "Setup")

        tk.Label(card, text="Function f(x)", bg=PANEL, fg=DIM,
                 font=("Segoe UI",9)).pack(anchor="w", padx=14, pady=10)
        self.v_func = tk.StringVar(value="exp(x)")
        self._ent(card, self.v_func).pack(fill="x", padx=14)

        row = tk.Frame(card, bg=PANEL)
        row.pack(fill="x", padx=14, pady=8)
        row.columnconfigure(0, weight=1); row.columnconfigure(1, weight=1)
        tk.Label(row, text="Lower (a)", bg=PANEL, fg=DIM,
                 font=("Segoe UI",9)).grid(row=0,column=0,sticky="w")
        tk.Label(row, text="Upper (b)", bg=PANEL, fg=DIM,
                 font=("Segoe UI",9)).grid(row=0,column=1,sticky="w",padx=6)
        self.v_a = tk.StringVar(value="0")
        self.v_b = tk.StringVar(value="2")
        self._ent(row, self.v_a).grid(row=1,column=0,sticky="ew")
        self._ent(row, self.v_b).grid(row=1,column=1,sticky="ew",padx=6)

        tk.Label(card, text="Intervals (n)", bg=PANEL, fg=DIM,
                 font=("Segoe UI",9)).pack(anchor="w", padx=14, pady=8)
        self.v_n = tk.StringVar(value="100")
        self._ent(card, self.v_n).pack(fill="x", padx=14)
        tk.Label(card, text="Higher values = more accuracy",
                 bg=PANEL, fg="#AAAACC", font=("Segoe UI",8)
                 ).pack(anchor="w", padx=14, pady=2)

        bf = tk.Frame(card, bg=PANEL)
        bf.pack(fill="x", padx=14, pady=12)
        tk.Button(bf, text="✦  Calculate",
                  bg=PURPLE, fg="white", font=("Segoe UI",11,"bold"),
                  relief="flat", cursor="hand2", pady=9,
                  activebackground="#4A4AC5",
                  command=self._calc).pack(fill="x")
        tk.Button(bf, text="Clear",
                  bg="#EBEBF5", fg=DIM, font=("Segoe UI",9),
                  relief="flat", cursor="hand2", pady=5,
                  command=self._clear).pack(fill="x", pady=5)

        self.v_err = tk.StringVar()
        tk.Label(card, textvariable=self.v_err, bg=PANEL, fg=RED,
                 font=("Segoe UI",9), wraplength=260, justify="left"
                 ).pack(anchor="w", padx=14, pady=6)

        # Quick examples
        card2 = self._card(parent)
        card2.pack(fill="x", pady=10, padx=2)
        self._hdr(card2, "Quick Examples")

        for label, expr, a, b, icon in [
            ("x²",     "x**2",   "0", "1",       "📐"),
            ("sin(x)", "sin(x)", "0", "3.14159", "~"),
            ("e^x",    "exp(x)", "0", "2",       "e"),
            ("1/x",    "1/x",    "1", "5",       "/"),
        ]:
            row2 = tk.Frame(card2, bg="#F5F5FC", cursor="hand2")
            row2.pack(fill="x", padx=8, pady=2)
            tk.Label(row2, text=icon, bg="#E8E8F8", fg=PURPLE,
                     font=("Segoe UI",11), width=3).pack(side="left", padx=6, pady=6)
            inner = tk.Frame(row2, bg="#F5F5FC")
            inner.pack(side="left", padx=4)
            tk.Label(inner, text=label, bg="#F5F5FC", fg=TEXT,
                     font=("Segoe UI",10,"bold")).pack(anchor="w")
            tk.Label(inner, text=f"[{a}, {b}]", bg="#F5F5FC", fg=DIM,
                     font=("Segoe UI",9)).pack(anchor="w")
            cb = lambda e, ex=expr, aa=a, bb=b: self._load(ex, aa, bb)
            for w in [row2, inner] + list(inner.winfo_children()):
                w.bind("<Button-1>", cb)

        tk.Frame(card2, height=4, bg=PANEL).pack()

        # Info
        info = tk.Frame(parent, bg="#EEF0FF", padx=10, pady=8)
        info.pack(fill="x", padx=2, pady=6)
        tk.Label(info, text="Supported Functions:", bg="#EEF0FF",
                 fg=PURPLE, font=("Segoe UI",9,"bold")).pack(anchor="w")
        for t in ["sin, cos, tan, exp, log, sqrt, abs",
                  "Constants: pi, e", "Use x as the variable"]:
            tk.Label(info, text=t, bg="#EEF0FF", fg="#4B4B80",
                     font=("Segoe UI",9)).pack(anchor="w")

    # ── Right panel ───────────────────────────────────────
    def _build_right(self, parent):
        # Ready placeholder
        self._ready = tk.Frame(parent, bg=PANEL)
        self._ready.pack(fill="both", expand=True)
        c = tk.Frame(self._ready, bg=PANEL)
        c.place(relx=.5, rely=.42, anchor="center")
        ib = tk.Frame(c, bg=BG, width=68, height=68)
        ib.pack(); ib.pack_propagate(False)
        tk.Label(ib, text="∫", bg=BG, fg=PURPLE,
                 font=("Segoe UI",30,"bold")).place(relx=.5,rely=.5,anchor="center")
        tk.Label(c, text="Ready to Calculate", bg=PANEL, fg=TEXT,
                 font=("Segoe UI",16,"bold")).pack(pady=12)
        tk.Label(c, text="Enter your function and bounds,\nthen click Calculate",
                 bg=PANEL, fg=DIM, font=("Segoe UI",10), justify="center").pack()

        # Results area (hidden until calculated)
        self._res_frame = tk.Frame(parent, bg=BG)
        self.v_title = tk.StringVar()
        tk.Label(self._res_frame, textvariable=self.v_title,
                 bg=BG, fg=TEXT, font=("Segoe UI",13,"bold")
                 ).pack(pady=10)
        self._cards_frame = tk.Frame(self._res_frame, bg=BG)
        self._cards_frame.pack(fill="both", expand=True)

    def _show_results(self, display):
        self._ready.pack_forget()
        self._res_frame.pack(fill="both", expand=True)
        self.v_title.set(f"Computing   ∫ {display} dx")

        for w in self._cards_frame.winfo_children():
            w.destroy()

        # 2x2 grid of result cards using simple Labels
        self._cards_frame.columnconfigure(0, weight=1)
        self._cards_frame.columnconfigure(1, weight=1)

        for i, (name, color, desc) in enumerate(METHODS):
            r, c = divmod(i, 2)

            val = self._results.get(name, float("nan"))
            val_s = f"{val:.6f}" if not math.isnan(val) else "Error"

            outer = tk.Frame(self._cards_frame, bg=color, highlightthickness=0)
            outer.grid(row=r, column=c, sticky="ew", padx=5, pady=5)

            tk.Label(outer, text=name, bg=color, fg="white",
                     font=("Segoe UI",11,"bold"),
                     anchor="w", padx=12, pady=8).pack(fill="x")
            tk.Label(outer, text=desc, bg=color, fg="#CCCCEE",
                     font=("Segoe UI",8),
                     anchor="w", padx=12, pady=6).pack(fill="x")

            body = tk.Frame(outer, bg=PANEL)
            body.pack(fill="both", expand=True)

            tk.Label(body, text=val_s,
                     bg=PANEL, fg=TEXT,
                     font=("Segoe UI",22,"bold"),
                     anchor="w", padx=12, pady=10
                     ).pack(fill="x")
            tk.Label(body, text="Result value",
                     bg=PANEL, fg=DIM,
                     font=("Segoe UI",8),
                     anchor="w", padx=12
                     ).pack(fill="x")

            btns = tk.Frame(body, bg=PANEL)
            btns.pack(anchor="w", padx=10, pady=8)

            tk.Button(btns, text="View Solution",
                      bg=PANEL, fg=color,
                      font=("Segoe UI",9,"bold"),
                      relief="solid", bd=1, cursor="hand2",
                      padx=8, pady=3,
                      command=lambda n=name, cl=color: self._trail(n, cl)
                      ).pack(side="left", padx=5)

            tk.Button(btns, text="Plot",
                      bg=PANEL, fg=DIM,
                      font=("Segoe UI",9),
                      relief="solid", bd=1, cursor="hand2",
                      padx=8, pady=3,
                      command=lambda n=name, cl=color: self._plot(n, cl)
                      ).pack(side="left")

    # ── Actions ───────────────────────────────────────────
    def _calc(self):
        self.v_err.set("")
        raw = self.v_func.get().strip()
        if not raw:
            self.v_err.set("Please enter a function."); return
        try:
            a = float(self.v_a.get())
            b = float(self.v_b.get())
            n = int(self.v_n.get())
        except ValueError:
            self.v_err.set("Invalid a, b, or n."); return
        if a >= b: self.v_err.set("a must be < b."); return
        if n < 2:  self.v_err.set("n must be >= 2."); return

        expr = clean(raw)
        try:
            make_f(expr)((a+b)/2)
        except Exception as ex:
            self.v_err.set(f"Function error: {ex}"); return

        self._computed = {"raw": raw, "expr": expr, "a": a, "b": b, "n": n}
        self._results  = {}
        f = make_f(expr)
        for name, fn in FUNCS.items():
            try:    self._results[name] = fn(f, a, b, n)
            except: self._results[name] = float("nan")

        self._show_results(f"{a} to {b}  {raw}")

    def _trail(self, method, color):
        c = self._computed
        res = self._results.get(method, 0.0)
        try:
            lines = TRAIL_BUILDERS[method](
                c["raw"], c["expr"], c["a"], c["b"], c["n"], res)
        except Exception:
            lines = [("b", traceback.format_exc())]
        TrailWindow(self, method, color, lines)

    def _plot(self, method, color):
        c = self._computed
        ChartWindow(self, c["raw"], c["expr"],
                    c["a"], c["b"], c["n"], method, color)

    def _clear(self):
        self.v_func.set(""); self.v_a.set("")
        self.v_b.set(""); self.v_n.set("100")
        self.v_err.set("")
        self._results = {}; self._computed = {}
        self._res_frame.pack_forget()
        self._ready.pack(fill="both", expand=True)

    def _load(self, expr, a, b):
        self.v_func.set(expr); self.v_a.set(a); self.v_b.set(b)

    # ── Helpers ───────────────────────────────────────────
    def _card(self, parent):
        return tk.Frame(parent, bg=PANEL,
                        highlightbackground=BORDER, highlightthickness=1)

    def _hdr(self, parent, text):
        tk.Label(parent, text=text, bg=PANEL, fg=TEXT,
                 font=("Segoe UI",11,"bold"),
                 anchor="w", padx=14, pady=9).pack(fill="x")
        tk.Frame(parent, bg=BORDER, height=1).pack(fill="x")

    def _ent(self, parent, var):
        wrap = tk.Frame(parent, bg="#F0F0FA",
                        highlightbackground=BORDER, highlightthickness=1)
        e = tk.Entry(wrap, textvariable=var,
                     bg="#F0F0FA", fg=TEXT,
                     font=("Segoe UI",11), relief="flat",
                     bd=0, insertbackground=PURPLE, highlightthickness=0)
        e.pack(fill="x", padx=6, pady=6)
        e.bind("<FocusIn>",  lambda _: wrap.config(highlightbackground=PURPLE))
        e.bind("<FocusOut>", lambda _: wrap.config(highlightbackground=BORDER))
        return wrap


if __name__ == "__main__":
    app = App()
    app.mainloop()