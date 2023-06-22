using LatticeModels
using JLD
using Ising2D
using Statistics

n_sweep = 2500           # number of sweeps between sampling
steps = 5000          # number of sweeps to thermalize
temps   = 2.7:-0.05:1.8  #temperatures to sample

function graf(;L)
    mt = Float64[] #Acá van a ir todos los puntos de la magnetización para cada T, los adjunto.
    xt = Float64[]
    IS = Ising(SQLattice_periodic,L,L)
    for T in temps      #loop en la temperatura, termaliza
        IS = Ising(SQLattice_periodic,L,L,T)
        Wolff!(IS,steps)
        m1=Float64[]
        for i=1:n_sweep
            Wolff!(IS,n_sweep)
            push!(m1,magnetization_ising2d(IS))                  
        end
        ma_ave = Statistics.mean(m1)  # Statistics.mean(m1)
        susceptibility = ((Lsize)^2)*(Statistics.var(m1)/T)
                
        push!(mt,ma_ave)
        push!(xt,susceptibility)
        
    end
     
    return collect(temps),xt
    
end



using JLD
using LatticeModels
using BioStatPhys
using DelimitedFiles

function compute_times(;L,ntimes,mlen)
    IS = Ising(SQLattice_periodic,L,L)
    conf = load("/Users/christian/SQconf_Tc_L$(L).jld","IS.σ")
    IS.σ .= conf
    set_energy_mag!(IS)
    set_temperature!(IS,Ising_SQ_critical_temperature)
    times=zeros(Float64,ntimes)
    fail = 0
    for i ∈ eachindex(times)
        try
            _,M = Metropolis!(IS,steps=mlen,save_interval=1)
            M = transpose(M)
            C=BioStatPhys.time_correlation_tw_direct(M,connected=true,i0=nothing,Xmean=zeros(size(M)),normalized=false)
            # C=time_correlation(M,connected=true,normalized=true,i0=1)
           
            times[i]=correlation_time_spectral(C,1)
        catch
            fail += 1
        end
    end
    if fail>0 @warn "Failed $(fail)" end

    writedlm("times_$(L)_$(ntimes)_$(mlen).txt",times)
    return times
end
