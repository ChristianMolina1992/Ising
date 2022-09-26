using Random: default_rng, seed!
using RandomNumbers
using StaticArrays
using Plots
using BenchmarkTools

const β_crit = log(1+sqrt(2))/2
rand_ising2d(m, n=m) = rand(Int8[-1, 1], m, n)

seed!(4649)
s₀ = rand_ising2d(100);


function energy_density_ising2d(s)
    m, n = size(s)
    E = 0.0
    @inbounds begin
        for j in 1:n, i in 1:m-1
            E -= s[i,j]*s[i+1,j]
        end
        for j in 1:n
            E -= s[m,j]*s[1,j]
        end
        for j in 1:n-1, i in 1:m
            E -= s[i,j]*s[i,j+1]
        end
        for i in 1:m
            E -= s[i,m]*s[i,1]
        end
    end
    return E/(m*n)
end

function ising2d_ifelse!(s, β, niters, rng=default_rng())
    m, n = size(s)
    prob = @SVector [exp(-2*β*k) for k in -4:4]
    @fastmath @inbounds @simd for iter in 1:niters
        for j in 1:n 
            for i in 1:m
                let NN = s[ifelse(i == 1, m, i-1), j],  #ifelse(condition::Bool, x, y)Return x if condition is true, otherwise return y
                    SS = s[ifelse(i == m, 1, i+1), j],
                    WW = s[i, ifelse(j == 1, n, j-1)],
                    EE = s[i, ifelse(j == n, 1, j+1)],
                    CT = s[i, j]
                    k = CT * (NN + SS + WW + EE)
                    s[i,j] = ifelse(rand(rng) < prob[k+5], -CT, CT)
                end
            end
        end
    end
end

s = copy(s₀)
rng = default_rng()
ising2d_ifelse!(s, β_crit, 10, rng)

s = copy(s₀)
seed!(4649)
rng = default_rng()
@time ising2d_ifelse!(s, β_crit, 10^5, rng)
heatmap(s; size=(200, 200), axis=false, colorbar=false)
function magnetization_ising2d(s)
    m, n = size(s)
    return sum(s)/(m*n)
end
s = rand_ising2d(300)
magnetization_ising2d(s)    # Acá me tira error in method definition: function Ising2D.energy_density_ising2d must be explicitly imported to be extended





using Printf
using Plots
using Ising2D

m       = 200   
n       = 200
niters  = 1000
n_sweep = 20             # number of sweeps between sampling
n_therm = 1000           # number of sweeps to thermalize
n_data  = 100            # number of data samples per temperature
temps   = 4.0:-0.3:0.1   # temperatures to sample

function sweep(p, T, s) # apply flip() to every site on the lattice
    for b = 1:p
        for i = 1:m
            for j = 1:n
                ising2d_ifelse!(s,T,niters, default_rng())
            end
        end
    end   
end




function grafico(temps)
    mt = Float64[] #Acá van a ir todos los puntos de la magnetización para cada T, los adjunto.
    et = Float64[]
    xt = Float64[]
    s = rand_ising2d(10)  # T=oo
    for T in temps
        sweep(n_therm, T, s)  #termaliza la red (loop sobre las temperaturas)

        m1 = Float64[]          #Acá irían todas las magnetizaciones diferentes para un T fijo
        e1 = Float64[]
        x1 = Float64[]
        for i in niters       #calcula la magnetizacion
            sweep(n_sweep, T, s)
            push!(m1,magnetization_ising2d(s))
            push!(e1,energy_density_ising2d(s))
        end
        ma_ave = sum(m1)  # Statistics.mean(m1)
        en_ave = sum(e1)  # Statistics.mean(e1)
        #susceptibility = Statistics.var(mt)/T
        push!(mt,ma_ave)
        push!(et,en_ave)
        # push!(x1,suscept)
        #println(ma_ave)
      
    end
    #plots(temps,ma_ave) # plot magnetization vs. temperature
    return temps,mt,et #susceptibility
end
