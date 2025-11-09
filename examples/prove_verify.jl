using BinaryFields, StatsBase
using Ligerito

# We provide multiple 'hardcoded' configurations such as
# hardcoded_config_20
# hardcoded_config_24
# hardcoded_config_28
# hardcoded_config_30
# Note that this _must_ match the verifier configuration and the size of the
# vector `poly` containing the polynomial coefficients.
config = Ligerito.hardcoded_config_24(BinaryElem32, BinaryElem128)
poly = rand(BinaryElem32, 2^24)

@info "Running with $(Threads.nthreads()) threads"
@time proof = prover(config, poly)

proof_size = sizeof(proof)
@info "Proof size: $(Base.format_bytes(proof_size))"

# We provide multiple 'hardcoded' configurations such as
# hardcoded_config_4_verifier()
# hardcoded_config_24_verifier()
# hardcoded_config_28_verifier()
# hardcoded_config_30_verifier()
# Note that this _must_ match the verifier configuration and the size of the
# vector `poly` containing the polynomial coefficients.
verifier_cfg = Ligerito.hardcoded_config_24_verifier()
@time verification_res = verifier(verifier_cfg, proof)
println("Verification result: ", verification_res)
@assert verification_res == true