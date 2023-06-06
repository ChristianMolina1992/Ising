using JLD
using LatticeModels
using BioStatPhys
using DelimitedFiles
using Plots 

function tauVStime(;L,nlong,ntimes,mlen)
    IS = Ising(SQLattice_periodic,L,L)
    conf = load("/Users/christian/SQconf_Tc_L$(L).jld","IS.σ")
    IS.σ = conf
    set_energy_mag!(IS)
    set_temperature!(IS,Ising_SQ_critical_temperature)
    epsilon=zeros(Float64,nlong)  
    times=zeros(Float64,ntimes)
    fail = 0
    
    i = 1
    pos = zeros(L^2, 2)
    for x = 1:L, y = 1:L
        pos[i, 1] = x
        pos[i, 2] = y
        i += 1
    end
     binning = distance_binning(pos, 1.0, rmin=0.0)
    
    for i ∈ eachindex(epsilon)   
    s=reshape(IS.σ,L^2)           
    r,C = space_correlation(binning,reshape(IS.σ,L^2),connected=true,normalized=true)
    epsilon[i] = correlation_length_r0(r,C)
    try
         _,M = Metropolis!(IS,steps=mlen,save_interval=1)
        M = transpose(M)
        D=BioStatPhys.time_correlation_tw_direct(M,connected=true,i0=1,Xmean=zeros(size(M)),normalized=true)
        times[i]=correlation_time_spectral(D,1)
        catch
          fail += 1
        end
    
    end
    
    if fail>0 @warn "Failed $(fail)" 
    end
    writedlm("tamaño_ro_$(L)_$(nlong)_$(mlen).txt",epsilon)
    
    writedlm("tamaño_tau_$(L)_$(ntimes)_$(mlen).txt",times)
    return epsilon,times
end



##p_100 = @time compute_times(L=100, ntimes=3000, mlen=550000)
