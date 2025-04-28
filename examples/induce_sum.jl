# create random multi-linear poly 
# encode it with ligero matrix 
# sample a few random challenges 
# partially evaluate that poly  
# sample a few random queries 
# induce a sumcheck given queries 
# make sure that the sumcheck is correct

using BinaryFields, MultilinearPoly, BinaryReedSolomon

@info "Running with $(Threads.nthreads()) threads"

# start with a random multi-linear polynomial
k = 14
poly = rand(BinaryElem16, 2^k)
f = MultiLinearPoly(poly)

# encode it with a ligero like cm
mat = encode_poly(poly)
rows = size(mat, 1)
cols = size(mat, 2)
cm = ligero_commit(poly, Int(rows / 4), cols, reed_solomon(BinaryElem16, Int(rows / 4), rows))
@assert mat == cm.mat
println("rows:", rows)
println("cols:", cols)

# prepare a basis evaluation 
n = Int(log2(cols))
log_basis_len = (k - n)
sks_vks = eval_sk_at_vks(2^log_basis_len, BinaryElem16)

# sample log(cols) random challenges
rs = rand(BinaryElem16, n)

f_eval = partial_eval(f, rs).evals 
println("len of f_eval: ", length(f_eval))


# now given the initial f and it's partial evaluation there should be a sumcheck relation holding 
# we can sample a few random queries for starting point 

# get random queries 
queries = sort(rand(1:rows, 148)) 
println("queries: ", queries)  

# get random separation challenge 
alpha = rand(BinaryElem16)

# open rows in queries: Notice that for query = 0 corresponds to a first row in the matrix
opened_rows = [vec(mat[q, :]) for q in queries]

# finally induce the sumcheck polynomial and enforced sum
basis_poly, enforced_sum = induce_sumcheck_poly(log_basis_len, sks_vks, opened_rows, rs, queries, alpha)

# evaluate f partially in those challenges, this is what hones prover commits to


println("len of basis_poly: ", length(basis_poly))

# now inner product of the partially evaluated f and the basis poly should be equal to the enforced sum 
inner_product = sum(f_eval .* basis_poly)
inner_product == enforced_sum

# v_ch = evaluate_lagrange_basis(rs)  
# println("len of v_ch: ", length(v_ch))



