function sanitize_transition(T::Matrix{<:Real}, I::Int64)
  @assert(#
    size(T,1) == size(T,2),
    "Transition probability matrix T must be square."
  )
  @assert(#
    I <= size(T,1),
    "Not a valid starting state."
  )
end

