using BinaryFields

config = hardcoded_config(BinaryElem32, BinaryElem128)
poly = rand(BinaryElem32, 2^24)

proof = prover(config, poly)

size_bytes = Base.summarysize(proof)
size_kb = size_bytes / 1024
println("Proof size: ", round(size_kb, digits=2), " KB")
# proof_size = sizeof(proof)
# @info "Proof size: $(Base.format_bytes(proof_size))"
