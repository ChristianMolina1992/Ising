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


function magnetization_ising2d(s)
    m, n = size(s)
    return sum(s)/(m*n)
end
s = rand_ising2d(300)
magnetization_ising2d(s)    # Acá me tira error in method definition: function Ising2D.energy_density_ising2d must be explicitly imported to be extended


function Var_de_m(s)
    m, n = size(s)
    return (std(magnetization_ising2d(s))^2)
end


function Var_de_e(s)
    m, n = size(s)
    return (std(energy_density_ising2d(s))^2)
end

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
                ising2d_ifelse!(s,T,niters, rng=default_rng())
            end
        end
    end   
end

function grafico()
    m1=(1:niters)        #Acá irían todas las magnetizaciones diferentes para un T fijo
    e1=(1:niters)
    X1=(1:niters)
    C1=(1:niters)
    mt = []              #Acá van a ir todos los puntos de la magnetización, los adjunto
    et = []
    Xt = []
    Ct = []
    s = rand_ising2d(200)
    for T in temps
        sweep(n_therm, T, s)  #termaliza la red
        m = Float[64]
        e = Float[64]
        X = Float[64]
        C = Float[64]
        push!(m1,magnetization_ising2d(s))
        push!(e1,energy_density_ising2d(s))
        push!(X1,Var_de_m(s))
        push!(C1,heat_capacity(e))        
        
        for i in niters       #calcula la magnetizacion
            sweep(n_sweep, T, s)
            magnetization_ising2d(s)
            energy_density_ising2d(s)
            Var_de_m(s)
            Var_de_e(s)
        end
        ma_ave = sum(m) 
        en_ave = sum(e)
        susceptibility = std(m)^2/temps
        heat_capacity = std(e)^2/temps
        push!(mt,ma_ave)
        push!(et,en_ave)
        @printf("%8.3f  %8.3f \n",ma_ave,en_ave,susceptibility,heat_capacity)
      
    end
    plots(temps,ma_ave) # plot magnetization vs. temperature
    plots(temps,ea_ave) # plot energy vs. temperature
    plots(temps,susceptibility)  #  susceptibility vs temperature
    plots(temps,heat_capacity)   # heat_capacity vs temperature
end


grafico()  #UndefVarError: ising2d_ifelse! not defined

