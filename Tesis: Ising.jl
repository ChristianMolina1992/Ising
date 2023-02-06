using Random: default_rng, seed!
using RandomNumbers
using StaticArrays
using Plots
using BenchmarkTools

const β_crit = log(1+sqrt(2))/2
rand_ising2d(m, n=m) = rand(Int8[-1, 1], m, n)
Lsize=100

seed!(4649)
s₀ = rand_ising2d(Lsize);


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


####PARA EL CALCULO DE CORRELACIONES#####


using Printf
using Plots
using Ising2D
using Statistics 
using DelimitedFiles


n_sweep = 1000            # number of sweeps between sampling
n_therm = 1500          # number of sweeps to thermalize
n_data  = 1000           # number of data samples per temperature
temps   = 2.4:-0.01:2.22   # temperatures to sample
Lsize= 10
function graf()
    mt = Float64[] #Acá van a ir todos los puntos de la magnetización para cada T, los adjunto.
    #et = Float64[]
    xt = Float64[]
    #ct = Float64[]
    s = rand_ising2d(Lsize) 
    for T in temps      #loop en la temperatura, termaliza
        ising2d_ifelse!(s,1/T,n_therm)
        m1=Float64[]
        e1=Float64[]
        writedlm("config_$(Lsize)_$(T).txt",s)
        for i=1:n_sweep
            ising2d_ifelse!(s,1/T,n_sweep)
            push!(m1,magnetization_ising2d(s))
            push!(e1,energy_density_ising2d(s))
                              
        end
        ma_ave = Statistics.mean(m1)  # Statistics.mean(m1)
        #en_ave = Statistics.mean(e1)  # Statistics.mean(e1)
        susceptibility = ((Lsize)^2)*(Statistics.var(m1)/T)
        #heat_capacity = Statistics.var(e1)/T^2
                
        push!(mt,ma_ave)
        #push!(et,en_ave)
        push!(xt,susceptibility)
        #push!(ct,heat_capacity)
        
    end
     
    return collect(temps),xt
    
end



using BioStatPhys
n_sweep = 10000
# Para usar esto: primero correr graf, luego llamar correlacion_temporal(T,Lsize) con una temperatura dada
@time function correlacion_temporal(T,L)
    spin = readdlm("config_$(L)_$(T).txt", Int8)
    m1=Float64[]
    for i=1:n_sweep
        ising2d_ifelse!(spin,1/T,1)
        push!(m1,magnetization_ising2d(spin))                                  
    end
    corr=time_correlation(m1,connected=true,normalized=true)
    time=collect(1:size(corr,1))
    #    return plot(time,corr,xlims=(10, 1000),xlabel = "Tiempo", ylabel = "C(t)",label="L=200, T = 3.0")
    return corr
end


# Para usar tcorr_vsT:
# cargar temps y Lsize
# luego llamar T,tau=tcorr_vsT()
# luego plotear tau vs T, tiene que dar un pico
function tcorr_vsT()
    tau=Float64[]
    for T in temps
        corr=correlacion_temporal(T,Lsize)
        push!(tau,correlation_time_spectral(corr,1))
    end
    return collect(temps),tau
end



#Para el calculo de la distribucion

n_sweep = 1000000
@time function distribuciones(T,Lsize)
    spin = readdlm("config_$(Lsize)_$(T).txt", Int8)  
    
    tiempo_relaj = Float64[]
    m1 = Float64[]
   
    for i=1:n_sweep
        ising2d_ifelse!(spin,1/T,1)
        push!(m1,magnetization_ising2d(spin))                                        
    end
    for t0 in 1:100
    corr = time_correlation(m1,connected=true)
    push!(tiempo_relaj,correlation_time_spectral(corr,t0*1000))  
    
    #En el código que me pasó, no funcionaba si ponía time_correlation(m1,connected=true,t0=t0, por lo que busque un por lo que busqué un poco más en la fuente de como funcionaba 
    y defenía a t0 = Deltat/1000, entonces por eso puse t0*1000 )
         
        
        end
        
        
     #Esta parte era para graficar el histograma, en principio el programa de arriba funciona correctamente, recién hoy lo logré hacer funcionar pero no sé si es correcto,
     esta sería la parte para graficar el histograma.
    struct Histogram
    his = Histogram(50, max = maximum(tiempo_relaj), min = minimun(tiempo_relaj))
    for tiempo in tiemporelaj 
    push!(his,tiempo) 
    tiempo,prob=prob(his)
    end
    end
    return plot(prob(his), tiempo)
    
        
end


#Con esto calculaba los histogramas

n_sweep = 200000
@time function distribuciones(T,Lsize)
    spin = readdlm("config_$(Lsize)_$(T).txt", Int8) 

    tiempo_relaj = Float64[]
    m1 = Float64[]

    for i=1:n_sweep
        ising2d_ifelse!(spin,1/T,1)
        push!(m1,magnetization_ising2d(spin))                                       
    end
    m1=transpose(m1)
    mmean=Statistics.mean(m1).*ones(size(m1))
    fail=0
    for t0 in 1:n_sweep÷5
        corr = BioStatPhys.time_correlation_tw_direct(m1,connected=true,i0=t0,Xmean=mmean,normalized=true)
        
        try
            tau=correlation_time_spectral(corr,1)
            push!(tiempo_relaj,tau)
        catch
            fail+=1
        end
    end
    if fail>0 @warn "Failed $(fail) tau" 
    end
    writedlm("tiempos_$(Lsize)_$(T).txt",tiempo_relaj)
    return collect(tiempo_relaj)
end
    


@time function histograma(T,Lsize)
time = readdlm("tiempos_$(T)_$(Lsize).txt",  '\t', '\n')
min,max=extrema(time)
    his = Histogram(100, max = max, min = min)
    for tiempo in time
        push!(his,tiempo)
    end
    return prob(his)
end

