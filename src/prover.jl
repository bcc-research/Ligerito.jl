using BinaryFields

struct LigeritoProver{T<:BinaryElem}
    dims::Vector{Tuple{Int, Int}}
end

function prover(prover::LigeritoProver{T}, poly::Vector{T}) where T <: BinaryElem
    # 0. commit to initial poly 
    # 1. get first k partial evaluations 
    # 2. commit to the claimed partial eval 
    # 3. receive S queries for the commitment
    # 4. induce a sumcheck based on queries
    # 5. instantiate a prover sumcheck instance  
    
    # recursive loop: 
        # for next k_i steps 
        # 0. receive a next evaluation point 
        # 1. run the folding step of the sumcheck 

        # after receiving last partial eval: 
        # 0. commit to next claimed partial evaluation 
        # 1. receive S queries for the commitment 
        # 2. induce a gluing sumcheck poly based on queries 
        # 3. 
end