function sanitize_transition(T::Matrix{<:Real}, I::Int64)
  @assert(#
    size(T,1) == size(T,2),
    "Transition probability matrix T must be square."
  )
  @assert(#
    I <= size(T,1),
    "Not a valid starting state."
  )
  @assert(#
    all([T[i,i] == 0 for i in 1:size(T,1)]),
    "All diagonal elements of the transition probability matrix T must be zero."
  )
end

