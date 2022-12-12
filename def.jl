mutable struct DR_data
    T :: Int64

    # market data
    D :: Array{Float64,1}
    h :: Array{Float64,1}
    s :: Array{Float64,1}

    # resource data
    cc :: Array{Float64,1}
    RD :: Array{Float64,1}
    RU :: Array{Float64,1}
    pmax :: Array{Float64,1}
    pmin :: Array{Float64,1}
    upmin :: Array{Float64,1}
    downmin :: Array{Float64,1}
    alpha :: Array{Float64,1}
    beta :: Array{Float64,1}
end
