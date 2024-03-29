using LatticeModels
using JLD
using Ising2D
using Statistics

nsample=1000
n_sweep = 250           # number of sweeps between sampling
steps = 5000          # number of sweeps to thermalize
temps   = 2.7:-0.05:1.8  #temperatures to sample

function graf(;L)
    mt = Float64[] #Acá van a ir todos los puntos de la magnetización para cada T, los adjunto.
    xt = Float64[]
    IS = Ising(SQLattice_periodic,L,L)
    for T in temps      #loop en la temperatura, termaliza
        set_temperature!(IS,T)
        Wolff!(IS,steps=steps)
        _,m1 = Wolff!(IS,steps=nsample*n_sweep,save_interval=n_sweep)
        ma_ave = Statistics.mean(m1)  # Statistics.mean(m1)
        susceptibility = (L^2)*(Statistics.var(m1)/T)
                
        push!(mt,ma_ave)
        push!(xt,susceptibility)
        
    end
     
    return collect(temps),xt
    
end



#ESTO SERIA LO NUEVO PARA CALCULAR CORRELACIONES TEMPORALES PROMEDIANDO SOBRE ORIGEN TEMPORAL

using JLD
using LatticeModels
using BioStatPhys
using DelimitedFiles

function compute_times_with(;L,ntimes,pasos)
    IS = Ising(SQLattice_periodic,L,L)
    conf = load("/Users/christian/SQconf_Tc_L$(L).jld","IS.σ")
    IS.σ .= conf
    set_energy_mag!(IS)
    set_temperature!(IS,Ising_SQ_critical_temperature)
    times=zeros(Float64,ntimes)
    fail = 0
    for t0 in pasos
        try
            _,M = Metropolis!(IS,steps=pasos,save_interval=1)
            C=time_correlation(M,connected=true)
                     
            times[i]=correlation_time_spectral(C,1)
        catch
            fail += 1
        end
    end
    if fail>0 @warn "Failed $(fail)" end

    #writedlm("times_$(L)_$(ntimes)_$(mlen).txt",times)
    return times
end
