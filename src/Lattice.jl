#=
In this file, I'll attempt to reproduce the code in LePage's Lattice QCD
for Novices paper in Julia. The code is meant to be reasonably faithful
to that in the paper apart from being in a different language.
=#

function update!(x::Vector{Float64}, ε::Float64)
    N = length(x)
	for j in 1:N
        old_x = x[j]
        old_Sj = S(j, x)
        x[j] = x[j] + uniform_float(-ε, ε)
        dS = S(j, x) - old_Sj
        if dS > 0 && exp(-dS) < uniform_float(0, 1)
            x[j] = old_x
        end
    end
end

function S(j::Int64, x::Vector{Float64}, a::Float64 = 0.5)
	N = length(x)
	jp = mod1(j+1, N)
	jm = mod1(j-1, N)
	# What is a? 
    return (a/2) * x[j]^2 + (x[j]/a) * (x[j] - x[jp] - x[jm])
end

function compute_G(x::Vector{Float64}, n::Int64)
    N = length(x)
    g = 0
    for j in 1:N
        g = g + x[j] * x[mod1(j+n, N)]
    end
    return g/N
end

function monte_carlo_average(x::Vector{Float64}, G, N_cor, ε, N_cf)
    N = length(x)
    for j in 1:N
        x[j] = 0
    end
    for j in 1:5*N_cor
        update!(x, ε)
    end
    for α in 1:N_cf
        for j in 1:N_cor
            update!(x, ε)
        end
        for n in 1:N
            G[α, n] = compute_G(x, n)
        end
    end
    for n in 1:N
        avg_G = 0
        for α in 1:N_cf
            avg_G += G[α, n]
        end
        avg_G = avg_G/N_cf
        print("G(", n, ") = ", avg_G, "\n")
    end
end

function bootstrap(G)
	N_cf = length(G)
    G_bootstrap = [G[Int64(floor(uniform_float(1, N_cf + 1)))] 
                   for i in 1:N_cf]
    return G_bootstrap
end

function bin(G, bin_size)
    G_binned = []
    for i in 1:bin_size:length(G)
        G_avg = 0
        for j in 1:bin_size
            G_avg += G[i+j-1]
        end
        append!(G_binned, G_avg / bin_size)
    end
    return G_binned
end

function uniform_float(lower, upper)
    return (upper - lower) * rand(Float64) + lower
end

function avg(G)
    return sum(G, dims=1) / length(G)
end

function sdev(G)
    abs(avg(G^2) - avg(G)^2)^(1/2)
end

function ΔE(G, a)
    avgG = avg(G)
    adE = log.(abs.(avgG[1:end-1] ./ avgG[2:end]))
    return adE/a
end

function main()
    N = 20
    N_cor = 20
    N_cf = 100
    a = 0.5
    ε = 1.4

    x = zeros(N)
    G = zeros(N_cf, N)

    monte_carlo_average(x, G, N_cor, ε, N_cf)

    print("Average G:\n", avg(G), "\n")
    print("Average G (binned):\n", avg(bin(G,4)), "\n")
    print("Average G (bootstrap):\n", avg(bootstrap(G)), "\n")

    print("ΔE:\n", ΔE(G, a))
    print(size(bootstrap(G)))
    print(typeof(ΔE(bootstrap(G), a)))
    print(size(ΔE(bootstrap(G), a)))
    print("ΔE (bootstrap):\n", ΔE(bootstrap(G), a))
end

main()
