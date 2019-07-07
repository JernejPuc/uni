# Satisfiability (SAT) reductions

The following implementations are provided:
 - reduction of n-queens to SAT
 - reduction of k-graph colouring to SAT

Once executed, problem constraints are produced in the [DIMACS](http://www.cs.ucsb.edu/~cappello/190B/lectures/sat/satformat.pdf) format. These can be used by [SAT solvers](https://msoos.github.io/cryptominisat_web/) to (eventually) find the appropriate solution.

## Usage

Both reductions are implemented in `Python3`.

### Examples

Each `.py` file contains a single function, which can be easily imported (assuming the file resides in your working directory):

```python
# n-queens
from nq_sat import reduce_nq_sat

n = 4
dimacs = reduce_nq_sat(4)
print(dimacs)
```

CLI functionality is provided for convenience. The following example automatically writes the constraints to `nq-4.txt`:

```sh
$ python nq_sat.py 4 nq-4.txt
```

For graph colouring, graph data (in DIMACS format) must first be read (in this case from `g1.txt`):

```sh
$ python kg_sat.py g1.txt 3 kg_g1-3.txt
```
