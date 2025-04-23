using BinaryFields, StatsBase

x = rand(BinaryElem32, 2^28)

@time xc = convert.(BinaryElem128, x)

function threaded_each(nt::Int, xs::AbstractVector, f::Function)
    n = length(xs)
    chunk_size = ceil(Int, n / nt)

    Threads.@sync for t in 1:nt
        Threads.@spawn begin
            start_idx = (t - 1) * chunk_size + 1
            end_idx = min(t * chunk_size, n)

            @inbounds for i in start_idx:end_idx
                f(xs[i])
            end
        end
    end
end


function convert_x(x::Vector{T}) where T <: BinaryElem
    converted = Vector{BinaryElem128}(undef, length(x))
    Threads.@threads for i in 1:length(x)
        converted[i] = convert(BinaryElem128, x[i])
    end
    return converted    
end

@time c = convert_x(x)

@time cc = threaded_each(Threads.nthreads(), x, x -> convert(BinaryElem128, x))

# @time begin
#     converted = Vector{BinaryElem128}(undef, length(x))
#     Threads.@threads for i in 1:length(x)
#         converted[i] = convert(BinaryElem128, x[i])
#     end
# end 

@show length(xc)