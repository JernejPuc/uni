# Evalution for the subset-sum problem

The performance of the following algorithms can be tested:
 - dynamic programming (unoptimised)
 - exhaustive search
 - greedy 2-approximation algorithm
 - FPTAS (fully-polynomial-time approximation scheme)

Note that while `dynamic` and `exhaustive` algorithms are exact, `greedy` and `fptas` are approximative.

## Usage

For performance reasons, the task is implemented in `Rust`. Algorithms are included in `subset_sum.rs`, which is to be used in `main.rs` like so:

```rust
mod subset_sum;
```

### Examples

Assuming you've created a `cargo` project, put `s1.txt` into the project directory and put both `main.rs` ans `subset_sum.rs` into its `src` subdirectory, you should be able to perform the following: 

```sh
$ cargo run --release dynamic s1.txt
    Finished release [optimized] target(s) in 0.09s
     Running `target\release\subset_sum.exe dynamic s1.txt`
dynamic:
n = 50
k = 1243228
r = 25372
t = 1.050792
```

```sh
$ cargo run --release exhaustive s1.txt
    Finished release [optimized] target(s) in 0.09s
     Running `target\release\subset_sum.exe exhaustive s1.txt`
exhaustive:
n = 50
k = 1243228
r = 25372
t = 0.012066
```

```sh
$ cargo run --release greedy s1.txt
    Finished release [optimized] target(s) in 0.10s
     Running `target\release\subset_sum.exe greedy s1.txt`
greedy:
n = 50
k = 1243228
r = 25372
t = 0.000044
```

For `fptas`, the `epsilon` value can be ommitted (defaults to `0.5`). If specified, it should be provided as the final argument:

```sh
$ cargo run --release fptas s1.txt 0.4
    Finished release [optimized] target(s) in 0.09s
     Running `target\release\subset_sum.exe fptas s1.txt 0.4`
fptas:
n = 50
k = 1243228
r = 25164
t = 0.000485
```

Optionally, you can write the results to a file, using:

```sh
$ cargo run --release dynamic s1.txt s1-dyn.txt
```

```sh
$ cargo run --release exhaustive s1.txt s1-exh.txt
```

```sh
$ cargo run --release greedy s1.txt s1-gre.txt
```

```sh
$ cargo run --release fptas s1.txt s1-fptas.txt 0.4
```
