using BinaryFields, StatsBase
using Ligerito

config = Ligerito.hardcoded_config_30(BinaryElem32, BinaryElem128)
poly = rand(BinaryElem32, 2^30)

@info "Running with $(Threads.nthreads()) threads"
@time proof = prover(config, poly)

proof_size = sizeof(proof)
@info "Proof size: $(Base.format_bytes(proof_size))"

verifier_cfg = Ligerito.hardcoded_config_30_verifier()
@time verification_res = verifier(verifier_cfg, proof)
println("Verification result: ", verification_res)
@assert verification_res == true