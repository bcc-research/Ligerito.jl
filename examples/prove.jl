using BinaryFields, StatsBase
using Ligerito

config = Ligerito.hardcoded_config_24(BinaryElem32, BinaryElem128)
poly = rand(BinaryElem32, 2^24)

@info "Running with $(Threads.nthreads()) threads"
@benchmark proof = prover(config, poly)

# proof_size = sizeof(proof)
# @info "Proof size: $(Base.format_bytes(proof_size))"

# @time verification_res = verifier(proof)
# println("Verification result: ", verification_res)
# @assert verification_res == true