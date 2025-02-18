# SAT9: A Survey Propagation SAT Solver

This code implements the Survey Propagation algorithm, a heuristic method for solving Boolean satisfiability (SAT) problems. It's based on message passing between clauses and literals, as described in the research literature on Belief Propagation and related techniques. This implementation is designed for educational purposes, to explore the behavior of Survey Propagation.

## What the Code Does

The program takes a SAT problem in conjunctive normal form (CNF) as input, attempts to find a satisfying assignment (or a good partial assignment), and outputs a reduced problem. The core algorithm iteratively updates two types of messages:

* **η (eta):** Clause-to-literal messages, representing the "urgency" for a literal to be true to satisfy the clause.
* **π (pi):** Literal-to-clause messages, representing the "flexibility" of a literal (how much it's *not* constrained by other clauses).

The algorithm also uses a "reinforcement" mechanism to encourage variables to commit to a particular value (true or false). After a number of iterations (or when convergence is reached), the solver identifies variables with a strong bias and fixes their values. It then simplifies the problem by removing satisfied clauses and literals made false by the fixed assignments. The result is a (hopefully smaller) residual SAT problem.

**Important Note:** This is a *heuristic* algorithm. It's not guaranteed to find a solution even if one exists. The output is a *partial* assignment and a *reduced* problem, not necessarily a complete solution. The reduced problem may still be unsatisfiable.

## Input Format

The input must be a text file in a simplified DIMACS CNF format:

* Each line represents a clause
* Literals within a clause are separated by spaces
* A literal is represented by a variable name (1 to 8 ASCII characters between `!` and `~`, inclusive, not beginning with `~`)
* A tilde (`~`) preceding a variable name indicates negation
* Lines starting with `~ ` are treated as comments and ignored
* Empty clauses (lines with no literals) are ignored but logged

**Example (`example.cnf`):**

```
~ x2 x3 ~x4
x1 x3 x4
~x1 x2 x4
~x1 ~x2 x3
~x2 ~x3 x4
~x1 ~x3 ~x4
x1 ~x2 ~x4
```

## Compilation

The code consists of three files: `sat9.c`, `gb_flip.c`, and `gb_flip.h`. `gb_flip.c` and `gb_flip.h` implement the GraphBase random number generator.

Compile using GCC:

```bash
gcc -c gb_flip.c -o gb_flip.o && gcc sat9.c gb_flip.o -o sat9
```

This creates an executable named `sat9`.

## Running the Solver

Run the solver with the input file redirected to standard input:

```bash
./sat9 < example.cnf
```

## Command-Line Options

The solver supports several command-line options to control its behavior:

| Option | Description | Default |
|--------|-------------|---------|
| `v<integer>` | Verbosity level. Higher values produce more output on stderr. | 1 |
| `h<integer>` | Hash table size (logarithm base 2). Adjusts the number of hash lists. | 8 |
| `b<integer>` | Input buffer size (in characters). Must be larger than the longest input line. | 1024 |
| `s<integer>` | Random seed. Use different seeds to get different runs. | 0 |
| `d<integer>` | Delta for periodic state reports (in mems). Reports progress every delta memory accesses. | 0 |
| `t<integer>` | Maximum number of iterations. | 1000 |
| `l<integer>` | Minimum number of iterations before reinforcement starts. | 5 |
| `c<integer>` | Confidence percentage (0-100). Variables with bias >= this are fixed. | 50 |
| `p<float>` | Damping factor (0.0 to 1.0) for reinforcement. | 0.99 |
| `e<float>` | Convergence threshold. Iterations stop when the maximum change in η is less than this. | 0.01 |

**Example with Options:**
```bash
./sat9 -v4 -s123 -t2000 -c80 < example.cnf
```

This sets verbosity to 4, the random seed to 123, the maximum iterations to 2000, and the confidence threshold to 80%.

## Output

The solver produces output on both standard output (stdout) and standard error (stderr).

### Standard Error (stderr):

* Input Summary: Reports the number of variables, clauses, and literals read
* Progress Reports: (If verbose is set appropriately or delta is used) Shows information about iterations, convergence, and memory usage
* Fixing Variables: Reports the number of variables fixed and the iteration number
* Unit Propagation: Reports if any additional variables were fixed by unit propagation
* Reduced Problem: Reports the number of clauses in the reduced problem and the name of the file where it's written
* Clause Size Distribution: Shows a histogram of the sizes of the clauses in the reduced problem
* Mems and Bytes: Total memory accesses (mems) and approximate memory usage (bytes)
* Errors: Reports any errors, such as file I/O problems or contradictions

### Standard Output (stdout):

* Partial Assignment: A space-separated list of literals that have been fixed. For example, `~x1 x3` means variable x1 is set to false and variable x3 is set to true
* `~~?`: Indicates that the solver found a contradiction (the problem is unsatisfiable, at least within the limits of the heuristic)

### Reduced Problem File:

The reduced problem is written to a file in `/tmp` with a name like `/tmp/sat9-<random_seed>.dat`. This file has the same format as the input file, representing the residual SAT problem after the partial assignment.

### Verbosity Levels (`v` option):

* 1 (Default): Basic information (input summary, fixing variables, reduced problem, mems/bytes)
* 2: Adds messages at the beginning of each iteration
* 4: Adds information about the maximum difference in η values during each iteration, and the reinforcement value
* 8: Adds information about resetting eta values (useful for debugging)
* 16: Prints a histogram of the π values after convergence or reaching the iteration limit
* 32: Prints all the π and η values for each variable after convergence or reaching the iteration limit, also printing the calculated probabilities for each literal to be 1, 0, or *

### Example with complete output

```bash
./sat9 -v32 < example.cnf > output.txt 2> error.txt
```

The command executes the sat9 program with verbosity level 32, taking input from example.cnf. Standard output is redirected to output.txt, and standard error is redirected to error.txt.

content of output.txt:
```
~x1 ~x2
```

content of error.txt:
```
(4 variables, 7 clauses, 21 literals successfully read)
beginning iteration 1
beginning iteration 2
beginning iteration 3
beginning iteration 4
beginning iteration 5
(fixing 2 variables after 5 iterations, e=0.00978189)
(unit propagation fixed 0 more variable(s))
Reduced problem of 2 clauses written on file /tmp/sat9-0.dat
 (1 3-clauses)
 (1 2-clauses)
Altogether 16356+38618 mems, 10768 bytes.
converged after 5 iterations.
variable      pi(v)        pi(~v)         1    0    *
      x1   0.001776(0)   0.997888(0)    0.00 0.99 0.00
      x2   0.002112(0)   0.997888(0)    0.00 0.99 0.00
      x3   0.668916(0)   0.355395(0)    0.42 0.23 0.35
      x4   0.488372(0)   0.511628(0)    0.26 0.25 0.49
```
