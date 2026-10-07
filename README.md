# Epidemic Spread Over N Weeks (C1)
**Computing Aᴺ with the Cayley–Hamilton Theorem to model a year of disease spread**

> Status: concept & design document. No code yet, the goal is to understand the problem first.

---

## 1. What the real-world problem means (simple terms)

Every person in a population is in one of three **states** in any given week:

| State | Meaning |
|---|---|
| **S** – Susceptible | healthy, can catch the disease |
| **I** – Infected | currently sick and infectious |
| **R** – Recovered | recovered, temporarily immune |

Each week, people **move between states with fixed probabilities**. For example, "30% of healthy people get infected this week" or "40% of infected people recover this week".

All those weekly probabilities are collected in one table, the **transition matrix A**.

**The questions we answer:**
1. After 1 week, 2 weeks, … 52 weeks, what fraction of people are S, I, R?
2. Can we get week 52 *directly*, without simulating 52 steps one by one?
3. In the long run, does the epidemic settle down? To what proportions? (**steady state**)

COVID modelling relies on this idea: SIR/SEIR-type models are repeated application of a "one-step" rule, i.e. matrix powers.

---

## 2. The exact mathematics

### 2.1 Markov chain form
Let the population state at week k be the column vector

```
x_k = [ S_k , I_k , R_k ]^T      (fractions, S_k + I_k + R_k = 1)
```

One week of evolution:

```
x_{k+1} = A · x_k      ⇒      x_n = Aⁿ · x_0
```

So everything reduces to computing **Aⁿ** (for n = 52 here).

### 2.2 Why not just multiply 52 times?
You *can*, and we will do it as a check. But Cayley–Hamilton gives a way to express **any** power of A using only I, A, A². That means a closed form in n, cheap for n = 52, 520 or 5,200,000.

### 2.3 Characteristic equation
For a 3×3 matrix A, the characteristic polynomial is

```
p(λ) = det(λI − A) = λ³ − τ λ² + σ λ − Δ
```

where
- τ = trace(A) = a₁₁ + a₂₂ + a₃₃
- σ = sum of the three principal 2×2 minors
- Δ = det(A)

### 2.4 Cayley–Hamilton theorem
> Every square matrix satisfies its own characteristic equation.

```
A³ − τ A² + σ A − Δ I = 0      ⇒      A³ = τ A² − σ A + Δ I
```

So A³ (and therefore every higher power) can be reduced to a combination of I, A, A².

### 2.5 General formula for Aⁿ
```
Aⁿ = c₀(n)·I + c₁(n)·A + c₂(n)·A²
```

The coefficients come from dividing λⁿ by p(λ): λⁿ = q(λ)·p(λ) + (c₀ + c₁λ + c₂λ²). Two ways to get them:

**Method A – recurrence (simplest to code).** Start with A⁰ = I, i.e. (c₀,c₁,c₂) = (1,0,0). Multiply by A and replace A³ using Cayley–Hamilton:

```
c₀' = Δ·c₂
c₁' = c₀ − σ·c₂
c₂' = c₁ + τ·c₂
```

Repeat n times. This only updates 3 numbers per step instead of multiplying 3×3 matrices.

**Method B – polynomial reduction by repeated squaring.** Compute λⁿ mod p(λ) with square-and-reduce in about log₂(n) steps (n = 52 needs ~6 squarings).

**Method C – eigenvalues (cross-check).** Solve p(λ) = 0 for the roots λ₁, λ₂, λ₃. Then c₀,c₁,c₂ solve
`λᵢⁿ = c₀ + c₁λᵢ + c₂λᵢ²` for i = 1,2,3 (a 3×3 linear system).

### 2.6 Important property of our matrix
A is **column-stochastic** (each column sums to 1, probabilities of leaving a state). Consequences:
- **λ = 1 is always an eigenvalue**, so p(λ) = (λ − 1)(λ² − (τ−1)λ + Δ). Only a quadratic remains.
- All other eigenvalues have |λ| ≤ 1, so Aⁿ converges to a fixed matrix (if |λ₂|, |λ₃| < 1).

### 2.7 Steady state
The steady state π is the vector that does not change:

```
A π = π      with      π_S + π_I + π_R = 1
```

It is the eigenvector for λ = 1, and also the limit of Aⁿx₀ as n → ∞ (every column of A^∞ equals π). The speed of convergence is set by the **second-largest eigenvalue modulus** max(|λ₂|, |λ₃|).

---

## 3. Input variables

| Symbol | Meaning | Typical value |
|---|---|---|
| **β** | weekly prob. a Susceptible becomes Infected (infection rate) | 0.30 |
| **γ** | weekly prob. an Infected becomes Recovered (recovery rate) | 0.40 |
| **δ** | weekly prob. a Recovered loses immunity → Susceptible (waning immunity) | 0.10 |
| **x₀** | initial state (S₀, I₀, R₀) | (0.99, 0.01, 0) |
| **n** | number of weeks | 52 |
| **N** (optional) | population size, to convert fractions to head-counts | 1,000,000 |

Constraints: 0 ≤ β, γ, δ ≤ 1 and S₀ + I₀ + R₀ = 1.

**Note on δ.** In the textbook SIR model, Recovered is a one-way trap (δ = 0), and the steady state is the trivial "everyone Recovered". Adding δ > 0 (an **SIRS** model, which matches COVID's waning immunity) gives a non-trivial steady state. We build the general model and treat δ = 0 as a special case.

---

## 4. How the model is formed

Columns = "from", rows = "to" (so x_{k+1} = A x_k):

```
            from S     from I     from R
to S   [   1 − β        0          δ    ]
to I   [    β         1 − γ        0    ]
to R   [    0           γ        1 − δ  ]
```

Reading it: a Susceptible stays S with 1−β or becomes I with β; an Infected stays I with 1−γ or becomes R with γ; a Recovered stays R with 1−δ or returns to S with δ. Every column sums to 1.

### Worked example (β = 0.3, γ = 0.4, δ = 0.1)

```
A = [ 0.7  0    0.1 ]
    [ 0.3  0.6  0   ]
    [ 0    0.4  0.9 ]
```

- τ = 0.7 + 0.6 + 0.9 = **2.2**
- σ = (0.7·0.6) + (0.7·0.9) + (0.6·0.9) = 0.42 + 0.63 + 0.54 = **1.59**
- Δ = det A = **0.39**

Characteristic equation:

```
λ³ − 2.2λ² + 1.59λ − 0.39 = 0
(λ − 1)(λ² − 1.2λ + 0.39) = 0
```

Eigenvalues: **λ₁ = 1**, **λ₂,₃ = 0.6 ± 0.1732 i** (complex pair, |λ| = √0.39 ≈ 0.6245).

What this tells us:
- Complex eigenvalues → the infected fraction shows **damped oscillation** (epidemic waves that fade).
- |λ| ≈ 0.62 < 1, and 0.6245⁵² ≈ 10⁻¹¹, so by week 52 the system has essentially reached the steady state.

Cayley–Hamilton reduction here: A³ = 2.2A² − 1.59A + 0.39I.

**First step by hand (check):** x₀ = (0.99, 0.01, 0) gives x₁ = (0.693, 0.303, 0.004). The entries sum to 1, as they should.

**Steady state:** flow balance βS = γI = δR gives π ∝ (1/β, 1/γ, 1/δ):

```
π = (1/0.3, 1/0.4, 1/0.1) / 15.833 ≈ (0.2105, 0.1579, 0.6316)
```

So long-run: ≈21% Susceptible, ≈16% Infected, ≈63% Recovered. Check: A·π = π ✓.

**Special case δ = 0 (plain SIR):** eigenvalues are 1, 1−β, 1−γ and the steady state is (0, 0, 1), i.e. everyone eventually recovers.

### Honest limitation (mention in your report)
In real epidemics β depends on how many people are infected (more infected → higher chance of catching it), which makes the model nonlinear. Our model uses a **constant** β, so it is a linear Markov chain and is exactly solvable with matrix powers. This is why it is a good first model and a building block for more realistic ones.

---

## 5. Algorithm to implement in Python

**Step 1: Build A** from β, γ, δ. Validate inputs (rates in [0,1], column sums = 1, x₀ sums to 1).

**Step 2: Characteristic equation.** Compute τ, σ, Δ. Get the polynomial coefficients [1, −τ, σ, −Δ].

**Step 3: Cayley–Hamilton power.** Compute c₀, c₁, c₂ for n (recurrence or repeated squaring). Form
`Aⁿ = c₀I + c₁A + c₂A²`.

**Step 4: Verify.** Compare with brute-force repeated multiplication and `numpy.linalg.matrix_power`. The difference should be ~10⁻¹⁵.

**Step 5: Evolve the epidemic.** For k = 0 … 52 compute x_k = A^k x₀ (using the Cayley–Hamilton coefficients for each k) and store the results in a table.

**Step 6: Eigen-analysis.** Find the roots of p(λ) (`numpy.roots`), check that λ = 1 is present, report |λ₂|, and estimate the convergence time.

**Step 7: Steady state.** Solve (A − I)π = 0 with Σπ = 1 (replace one equation by the sum constraint), and compare with the closed form (1/β, 1/γ, 1/δ) normalized, and with A^52 x₀.

**Step 8: Report and plot.** Show the week-by-week table, the S/I/R curves, the peak infection week, and the convergence to π.

**Complexity talking point:** brute force costs n matrix multiplications (O(n·3³)). The Cayley–Hamilton recurrence costs O(n) scalar updates on 3 numbers, and with repeated squaring O(log n).

---

## 6. What the final demo application should do

**Inputs (user-adjustable):** β, γ, δ, initial infected fraction, population N, number of weeks.

**Outputs:**
1. The matrix A, its characteristic polynomial, and the eigenvalues (with a note on whether they are complex → oscillation).
2. The Cayley–Hamilton expansion: A³ = τA² − σA + ΔI, then **A⁵²** as a 3×3 matrix, with c₀, c₁, c₂ displayed.
3. **Verification panel:** Cayley–Hamilton result vs. brute force vs. NumPy, with the error.
4. A **52-week table** of S, I, R (fractions and head-counts).
5. **Plots:** S/I/R vs. week (with the steady-state values as dashed lines) and optionally a simplex/triangle trajectory.
6. **Key facts:** peak infection week and size, week when the system is within 1% of steady state.
7. **Steady state π**, with a check A·π = π.
8. **What-if comparison:** e.g. δ = 0 (plain SIR) versus δ > 0 (SIRS), and a higher β scenario (a "lockdown vs. no lockdown" story).

**Suggested interface:** a command-line menu first (works anywhere), then an optional Streamlit page with sliders for β, γ, δ.

---

## 7. Suggested project structure

```
epidemic_matrix/
├── README.md
├── model.py          # build A, validation
├── cayley_hamilton.py# τ, σ, Δ, coefficients, A^n
├── analysis.py       # eigenvalues, steady state, peak week
├── demo.py           # CLI / Streamlit app and plots
└── tests/            # compare against numpy.linalg.matrix_power
```

## 8. Report checklist
- [ ] Define states, matrix A and meaning of β, γ, δ
- [ ] Derive the characteristic equation and show λ = 1 is a root
- [ ] State and apply Cayley–Hamilton; derive the c₀,c₁,c₂ recurrence
- [ ] Compute A⁵² and verify numerically
- [ ] Show the 52-week evolution
- [ ] Derive the steady state both algebraically and as limₙ Aⁿ
- [ ] Discuss complex eigenvalues (waves) and limitations (constant β)
