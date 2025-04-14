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

export setup 