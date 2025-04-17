using BinaryFields, StatsBase

config = hardcoded_config(BinaryElem32, BinaryElem128)
poly = rand(BinaryElem32, 2^24)

@info "Running with $(Threads.nthreads()) threads"
@time proof = prover(config, poly)

proof_size = sizeof(proof)
@info "Proof size: $(Base.format_bytes(proof_size))"