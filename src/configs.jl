# This can be done programmatically, but for now we just hardcode configs

function hardcoded_config_20(::Type{T}, ::Type{U}) where {T, U <: BinaryElem}
    recursive_steps = 1
    inv_rate = 4

    dims = Vector{Tuple{Int, Int}}(undef, recursive_steps)
    initial_dims = (2^14, 2^6)
    dims[1] = (2^10, 2^4)

    ks = Vector{Int}(undef, recursive_steps)
    initial_k = 6 
    ks[1] = 4 

    initial_reed_solomon = reed_solomon(T, initial_dims[1], initial_dims[1] * inv_rate)
    reed_solomon_codes = Vector{BinaryReedSolomon.ReedSolomonEncoding{U}}(undef, recursive_steps)
    for i in 1:recursive_steps
        reed_solomon_codes[i] = reed_solomon(U, dims[i][1], dims[i][1] * inv_rate)
    end

    return ProverConfig(recursive_steps, initial_dims, dims, initial_k, ks, initial_reed_solomon, reed_solomon_codes)
end

function hardcoded_config_20_verifier()
    recursive_steps = 1

    log_rs_dims = Vector{Int}(undef, recursive_steps)
    initial_dim = 14
    log_rs_dims[1] = 10

    ks = Vector{Int}(undef, recursive_steps)
    initial_k = 6 
    ks[1] = 4 

    return VerifierConfig(recursive_steps, initial_dim, log_rs_dims, initial_k, ks)
end

function hardcoded_config_24(::Type{T}, ::Type{U}) where {T, U <: BinaryElem}
    recursive_steps = 2
    inv_rate = 4

    dims = Vector{Tuple{Int, Int}}(undef, recursive_steps)
    initial_dims = (2^18, 2^6)
    dims[1] = (2^14, 2^4)
    dims[2] = (2^10, 2^4)

    ks = Vector{Int}(undef, recursive_steps)
    initial_k = 6 
    ks[1] = 4 
    ks[2] = 4 

    initial_reed_solomon = reed_solomon(T, initial_dims[1], initial_dims[1] * inv_rate)
    reed_solomon_codes = Vector{BinaryReedSolomon.ReedSolomonEncoding{U}}(undef, recursive_steps)
    for i in 1:recursive_steps
        reed_solomon_codes[i] = reed_solomon(U, dims[i][1], dims[i][1] * inv_rate)
    end

    return ProverConfig(recursive_steps, initial_dims, dims, initial_k, ks, initial_reed_solomon, reed_solomon_codes)
end

function hardcoded_config_24_verifier()
    recursive_steps = 2

    log_rs_dims = Vector{Int}(undef, recursive_steps)
    initial_dim = 18
    log_rs_dims[1] = 14
    log_rs_dims[2] = 10

    ks = Vector{Int}(undef, recursive_steps)
    initial_k = 6 
    ks[1] = 4 
    ks[2] = 4 

    return VerifierConfig(recursive_steps, initial_dim, log_rs_dims, initial_k, ks)
end

function hardcoded_config_28(::Type{T}, ::Type{U}) where {T, U <: BinaryElem}
    recursive_steps = 4 
    inv_rate = 4

    dims = Vector{Tuple{Int, Int}}(undef, recursive_steps)
    initial_dims = (2^22, 2^6)
    dims[1] = (2^19, 2^3)
    dims[2] = (2^16, 2^3)
    dims[3] = (2^13, 2^3)
    dims[4] = (2^10, 2^3)
 
    ks = Vector{Int}(undef, recursive_steps)
    initial_k = 6 
    ks[1] = 3 
    ks[2] = 3 
    ks[3] = 3 
    ks[4] = 3

    initial_reed_solomon = reed_solomon(T, initial_dims[1], initial_dims[1] * inv_rate)
    reed_solomon_codes = Vector{BinaryReedSolomon.ReedSolomonEncoding{U}}(undef, recursive_steps)
    for i in 1:4
        reed_solomon_codes[i] = reed_solomon(U, dims[i][1], dims[i][1] * inv_rate)
    end

    return ProverConfig(recursive_steps, initial_dims, dims, initial_k, ks, initial_reed_solomon, reed_solomon_codes)
end

function hardcoded_config_28_verifier()
    recursive_steps = 4

    log_rs_dims = Vector{Int}(undef, recursive_steps)
    initial_dim = 22
    log_rs_dims[1] = 19
    log_rs_dims[2] = 16
    log_rs_dims[3] = 13
    log_rs_dims[4] = 10

    ks = Vector{Int}(undef, recursive_steps)
    initial_k = 6 
    ks[1] = 3 
    ks[2] = 3 
    ks[3] = 3 
    ks[4] = 3

    return VerifierConfig(recursive_steps, initial_dim, log_rs_dims, initial_k, ks)
end

function hardcoded_config_30(::Type{T}, ::Type{U}) where {T, U <: BinaryElem}
    recursive_steps = 3 
    inv_rate = 4

    dims = Vector{Tuple{Int, Int}}(undef, recursive_steps)
    initial_dims = (2^23, 2^7)
    dims[1] = (2^19, 2^4)
    dims[2] = (2^15, 2^4)
    dims[3] = (2^11, 2^4)

    ks = Vector{Int}(undef, recursive_steps)
    initial_k = 7
    ks[1] = 4
    ks[2] = 4 
    ks[3] = 4 

    initial_reed_solomon = reed_solomon(T, initial_dims[1], initial_dims[1] * inv_rate)
    reed_solomon_codes = Vector{BinaryReedSolomon.ReedSolomonEncoding{U}}(undef, recursive_steps)
    for i in 1:recursive_steps
        reed_solomon_codes[i] = reed_solomon(U, dims[i][1], dims[i][1] * inv_rate)
    end

    return ProverConfig(recursive_steps, initial_dims, dims, initial_k, ks, initial_reed_solomon, reed_solomon_codes)
end

function hardcoded_config_30_verifier()
    recursive_steps = 3

    log_rs_dims = Vector{Int}(undef, recursive_steps)
    initial_dim = 23
    log_rs_dims[1] = 19
    log_rs_dims[2] = 15
    log_rs_dims[3] = 11

    ks = Vector{Int}(undef, recursive_steps)
    initial_k = 7
    ks[1] = 4
    ks[2] = 4 
    ks[3] = 4 

    return VerifierConfig(recursive_steps, initial_dim, log_rs_dims, initial_k, ks)
end