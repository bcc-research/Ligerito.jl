# Ligerito
This is an experimental implementation of Ligerito polynomial commitment scheme
described in [this paper](angeris.github.io/papers/ligerito.pdf) by Andrija
Novakovic and Guillermo Angeris.

**The code is not audited and is meant for research purposes only.**

# How to run this code
1. Clone the repository
2. Launch a Julia REPL in the cloned directory
3. Activate the environment (e.g., `using Pkg; Pkg.activate(".")`)
4. Instantiate all dependencies (via `Pkg.instantiate()`)
5. Run the example: (via `include("./examples/prove_verify.jl")`)

We provide multiple possible configurations. For more, see the
`examples/prove.jl` and `src/configs.jl` files.

# Note
We plan to upstream this into the `bcc-research/CryptoUtilitiesExamples`
repository and make this slightly more modular for use in future
implementations, but the current code serves as a reference and mirrors the
description of the algorithm in the Ligerito paper.