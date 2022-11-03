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
using Statistics 

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




using Printf
using Plots
using Ising2D
using Statistics 


n_sweep = 100             # number of sweeps between sampling
n_therm = 1000           # number of sweeps to thermalize
n_data  = 100            # number of data samples per temperature
temps   = 4.0:-0.3:0.1   # temperatures to sample

using Printf
using Plots
using Ising2D
using Statistics 


n_sweep = 100             # number of sweeps between sampling
n_therm = 1000           # number of sweeps to thermalize
n_data  = 100            # number of data samples per temperature
temps   = 4.0:-0.3:0.1   # temperatures to sample

function grafico()
    mt = Float64[] #Acá van a ir todos los puntos de la magnetización para cada T, los adjunto.
    et = Float64[]
    xt = Float64[]
    ct = Float64[]
    s = rand_ising2d(100) 
    for T in temps      #loop en la temperatura, termaliza
        ising2d_ifelse!(s,1/T,n_therm)
        m1=Float64[]
        e1=Float64[]
        
        for i=1:n_sweep
            ising2d_ifelse!(s,1/T,n_sweep)
            push!(m1,magnetization_ising2d(s))
            push!(e1,energy_density_ising2d(s))
                              
        end
        ma_ave = Statistics.mean(m1)  # Statistics.mean(m1)
        en_ave = Statistics.mean(e1)  # Statistics.mean(e1)
        susceptibility = Statistics.var(m1)/T
        heat_capacity = Statistics.var(e1)/T^2
                
        push!(mt,ma_ave)
        push!(et,en_ave)
        push!(xt,susceptibility)
        push!(ct,heat_capacity)
        
    end
     
    return collect(temps),xt
    
end
    

x=[4.0, 3.7, 3.4, 3.1, 2.8, 2.5, 2.2, 1.9, 1.6, 1.3, 1.0, 0.7, 0.4, 0.1]    #(m vs T)
y=[0.0026759999999999996, 0.00313, 0.002428, -0.0027800000000000004, -0.0038280000000000002, -0.006083999999999999, -0.78014, -0.9387600000000001, -0.9793999999999999, -0.994926, -0.9992779999999999, -0.9999819999999999, -1.0, -1.0]
plot(x, y, xlabel = "T" , ylabel = "Magnetization", seriestype = :scatter)

    
x=[4.0, 3.7, 3.4, 3.1, 2.8, 2.5, 2.2, 1.9, 1.6, 1.3, 1.0, 0.7, 0.4, 0.1]
y=[-0.5606319999999999, -0.6152639999999999, -0.6871439999999999, -0.7797919999999999, -0.9076639999999999, -1.108528, -1.5427079999999997, -1.8119079999999999, -1.9269800000000001, -1.9802840000000002, -1.997084, -1.9999440000000002, -2.0, -2.0]
plot(x, y, xlabel = "T", ylabel = "Energy", seriestype = :scatter)
    
