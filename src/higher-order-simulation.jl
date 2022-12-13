#using Gillespie
using Plots
using Random
using Distributions
using StaticArrays

"A type storing the status at the end of a call to `ssa`:
- **termination_status** : whether the simulation stops at the final time (`finaltime`) or early due to zero propensity function (`zeroprop`)
- **nsteps** : the number of steps taken during the simulation.
"
struct SSAStats
    termination_status::String
    nsteps::Int64
end

"A type storing the call to `ssa`:
- **x0** : a `Vector` of `Int64`, representing the initial states of the system.
- **F** : a `Function` or a callable type, which itself takes two arguments; x, a `Vector` of `Int64` representing the states, and parms, a `Vector` of `Float64` representing the parameters of the system.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`).
- **alg** : the algorithm used (`Symbol`, either `:gillespie`, `jensen`, or `tjc`).
- **tvc** : whether rates are time varying.
"
struct SSAArgs{X,Ftype,N,P}
    x0::X
    F::Ftype
    nu::N
    parms::P
    tf::Float64
    alg::Symbol
    tvc::Bool
end

"
This type stores the output of `ssa`, and comprises of:
- **time** : a `Vector` of `Float64`, containing the times of simulated events.
- **data** : a `Matrix` of `Int64`, containing the simulated states.
- **stats** : an instance of `SSAStats`.
- **args** : arguments passed to `ssa`.
"
struct SSAResult
    time::Vector{Float64}
    data::Matrix{Int64}
    stats::SSAStats
    args::SSAArgs
end


"
This function is a substitute for `StatsBase.sample(wv::WeightVec)`, which avoids recomputing the sum and size of the weight vector, as well as a type conversion of the propensity vector. It takes the following arguments:
- **w** : an `Array{Float64,1}`, representing propensity function weights.
- **s** : the sum of `w`.
- **n** : the length of `w`.
"
function pfsample(w::AbstractArray{Float64,1},s::Float64,n::Int64)
  t = rand() * s
  i = 1
  cw = w[1]
  while cw < t && i < n
    i += 1
    cw += w[i]
  end
  return i
end


"
Performs stochastic particle simulation algorithm.  
- **S** : m by n stoichiometry matrix with m metabolites and n reactions.
- **k** : n-length vector of kinetic parameters.
- **I** : Initial particle state.
- **tf** : the final simulation time.
"
function spsa(S::Matrix{<:Real},k::Vector{<:Real},I::Int64I::Int64,tf::<:Real)

  # Error checking
  @assert(all(k .> 0), "Kinetic parameter must be greater than zero.")
  @assert(size(S,2) == length(k), "Columns in S must equal length of k.")
  @assert(I ∈ 1:size(S,1), "Initial state I cannot exceed number of rows in S.")
  @assert(tf > 0, "Final simulation time must be greater than zero.")

  # Arguments
  #args = SSAArgs(x0,F,nu,parms,tf,:gillespie,false)

  # Total time array
  tt = Vector{Vector{Float64}}()

  # Inner time array
  ta = Vector{Float64}()
  push!(ta, 0.0)

  # Total particle trajectory array
  xt = Vector{Vector{Int64}}()

  # Inner particle trajectory array
  xa = Vector{Float64}()
  push!(xa, I)

  # Set up metabolite deficit array
  da = zeros(Int64, size(S,1))

  # Recursively simulate trajectories
  function spsa_inner(t, S, k, x)
    if t <= tf
      # Kinetic parameters and product states
      ks, ss, ps = F(S, k, Int64(x[end]))
      # Update time
      sumk = sum(ks)
      dt = rand(Exponential(1/sumk))
      push!(ta, ta[end] + dt)

      # Update trajectories
      j = pfsample(ks,sumk,length(ks))

      # Delete trajectory if current state matches a deficit token; remove token


      # Add token if reaction involves two or more substrates

      xx = [x for i in 1:length(j)]
      [push!(xx[i], ps[j][i]) for i in 1:length(ps[j])]

      # Continue simulating states
      for k ∈ ps[j]
        spsa_inner()
      end
    end
  end


    #if x isa SVector
      #@inbounds x[1] += nu[ev,:]
    #else
      #deltax = view(nu,ev,:)
      #for i in 1:nstates
        #@inbounds x[1,i] += deltax[i]
      #end
    #end
    #for xx in x
        #push!(xa,xx)
    #end


  #stats = SSAStats(termination_status,nsteps)
  #xar = transpose(reshape(xa,length(x),nsteps+1))
  #return SSAResult(ta,xar,stats,args)
end


# m[1]: Pseudo start
# m[2]: Glc
# m[3]: Pep
# m[4]: G6PHackathons sound like a scary  thing given my nervousness around meeting new people and dealing with very very new
# m[5]: Pyr
# m[6]: FGP
# m[7]: NADP
# m[8]: H2O
# m[9]: 6PGL
# m[10]: H
# m[11]: NADPH
# m[12]: 6PGC
# m[13]: CO2Hackathons sound like a scary  thing given my nervousness around meeting new people and dealing with very very new
# m[14]: Ru5P
# m[15]: Pseudo end

S = [#type
  -1 -1 -1 -1  0  0  0  0  0  0  0  0  0  0  0  0  1
   1  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0
   0  1  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  1  0 -1  1  0  0  0  0  0  0  0  0  0
   0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  1 -1 -1  1  0  0  0  0  0  0  0
   0  0  1  0  0  0  0  0 -1  1  0  0 -1  0  0  0  0
   0  0  0  1  0  0  0  0  0  0  0 -1  0  0  0  0  0
   0  0  0  0  0  0  0  0  1 -1  0 -1  0  0  0  0  0
   0  0  0  0  0  0  0  0  1 -1 -1  1  0  0  0  0  0
   0  0  0  0  0  0  0  0  1 -1  0  0  1 -1  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  1  0 -1  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  1  0  0 -1  0
   0  0  0  0  0  1  0  0  0  0  1  0  0  1  1  1 -1
]
k = [0.1, 0.4, 0.2, 0.3, 0.7, 0.4, 0.6, 0.4, 0.4, 0.45, 0.6, 0.1, 0.01, 0.1, 0.1, 0.9, 6]

function F(S::Matrix{<:Real}, k::Vector{<:Real}, i::Int64)
  cols = findall(S[i,:] .< 0)
  subs = [findall(S[:,j] .< 0) for j in cols]
  [subs[k] = subs[i][subs[i] .!= i] for k in 1:length(subs)]
  return k[cols], subs, [findall(S[:,j] .> 0) for j in cols]
end


Random.seed!(1234)
#result = gillespie(x0,F,nu,parms,tf)
#plot(result.data)






using Random


function simulate_stochastic_kinetics(#
  S::Matrix{<:Real},
  k::Vector{<:Real},
  I::Int64,
  t::Real,
  rng::AbstractRNG,
  c::Vector{<:Real}=ones(size(S,1))
)
  # More thorough error checking later
  @assert(size(S,2) == length(k))
  @assert(I ∈ 1:size(S,1))
  #check for closed-loop network

  X = Vector{Int64}()
  τ = Vector{Float64}()


  return 1
end



function ctmc(
    Q::Matrix{Float64},
    I::Int64, total_time::Float64, rng::AbstractRNG)
  @assert size(Q)[1] == size(Q)[2] # transition matrix is square
  @assert I ∈ 1:size(Q)[1]         # valid initial state

  # Total number of states
  N = size(Q)[1]

  # Check rows of Q sum to zero
  for i in 1:size(Q)[1]
    @assert sum(Q[i,:]) == 0
  end

  # Clear negative diagonal elements from infinitesimal generator matrix
  Q = copy(Q)
  for i in 1:size(Q)[1]
    Q[i,i] = 0
  end

  # Convert infinitesimal generator to probability matrix
  total_rates = sum(Q, dims = 2)
  P = Q ./ total_rates
  
  # Treat row vectors as categorical distributions
  dists = [Categorical(P[i, :]) for i in 1:N]

  # Initialize simulation and set the initial state
  X = Array{Int64,1}()
  τ = Array{Float64,1}()
  push!(X, I)
  push!(τ, 0.0)

  # Simulate embedded Markov chain
  while τ[end] < total_time
    push!(τ, τ[end] + rand(rng, Exponential(1 / total_rates[X[end]])))
    push!(X, rand(rng, dists[X[end]]))
  end

  # Return time, sequence of states, and probability transition matrix
  return τ, X, P
end

