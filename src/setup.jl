is_pow_two(x) = x & (x - 1) == 0

struct LigeroProofProperties
    log2_prob::Int
    inv_rate::Int
    base_field_size::Int
    ext_field_size::Int
end

n_queries(; log2_prob=100, inv_rate=4) = ceil(Int, -log2_prob/log2((1+1/inv_rate)/2))

# Assumes `N` is a power of 2
function opt_dims(N, prop::LigeroProofProperties)
    S = n_queries(; prop.log2_prob, inv_rate=prop.inv_rate)
    x = sqrt(prop.ext_field_size * N / (prop.base_field_size * S))
    n = 2 ^ ceil(Int, log2(x))
    return div(N, n), n
end

# basic setup
function setup(N::Int, base_field_size::Int, stop_crit::Int)
    initial_props = LigeroProofProperties(100, 4, base_field_size, 128)
    dims = []
    m, n = opt_dims(N, initial_props)
    push!(dims, (m, n))

    # after the first recursive call everything is in the ext field
    props = LigeroProofProperties(100, 4, 128, 128)
    while m > stop_crit
        m, n = opt_dims(m, props)
        push!(dims, (m, n))
    end

    return dims
end

function hardcoded_config(::Type{T}, ::Type{U}) where {T, U <: BinaryElem}
    # hardcode optimal dimensions for N = 2^24
    dims = Vector{Tuple{Int, Int}}(undef, 2)
    initial_dims = (2^18, 2^6)
    dims[1] = (2^14, 2^4)
    dims[2] = (2^10, 2^4)

    # hardcode amount of random variables we do per recursive step 
    ks = Vector{Int}(undef, 2)
    initial_k = 6 
    ks[1] = 4 
    ks[2] = 4 

    # prepare reed solomon codes 
    # 4 is hardcoded inv_rate for now
    initial_reed_solomon = reed_solomon(T, initial_dims[1], initial_dims[1] * 4)
    reed_solomon_codes = Vector{BinaryReedSolomon.ReedSolomonEncoding{U}}(undef, 2)
    for i in 1:2
        reed_solomon_codes[i] = reed_solomon(U, dims[i][1], dims[i][1] * 4)
    end

    recursive_steps = 2

    return ProverConfig(recursive_steps, initial_dims, dims, initial_k, ks, initial_reed_solomon, reed_solomon_codes)
end

export hardcoded_config