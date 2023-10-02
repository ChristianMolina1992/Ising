using JLD
using LatticeModels
using BioStatPhys
using DelimitedFiles

function prueba(;L, ntimes, mlen)
    IS = Ising(SQLattice_periodic, L, L;)
    conf = load("/Users/christian/SQconf_Tc_L$(L).jld", "IS.σ")
    IS.σ .= conf
    set_energy_mag!(IS)
    set_temperature!(IS, Ising_SQ_critical_temperature)
    times = zeros(Float64, ntimes)
    correlacion = zeros(Float64, ntimes)
    fail = 0
    for i ∈ eachindex(times)
        try
            _, M = Metropolis!(IS, steps = mlen, save_interval = 1)
            M = transpose(M)
            C= BioStatPhys.time_correlation_tw_direct(M, connected = true, i0 = 1, Xmean = zeros(size(M)), normalized = true)
            #println("Iteración $i: C = $C")
            times[i]=correlation_time_spectral(C,1)
            push!(correlacion,C)
        catch
            fail += 1
        end
    end

    if fail > 0
        @warn "Failed $(fail)"
    end

    #writedlm("times_$(L)_$(ntimes)_$(mlen).txt", times)
    return correlacion, times
end



using JLD
using LatticeModels
using BioStatPhys
using DelimitedFiles
using Plots 

function prueba_corr_rvs_C(;L,nlong,mlen)
    IS = Ising(SQLattice_periodic,L,L)
    conf = load("/Users/christian/SQconf_Tc_L$(L).jld","IS.σ")
    IS.σ = conf
    #matriz=zeros((L^2),longitud)
    set_energy_mag!(IS)
    set_temperature!(IS,Ising_SQ_critical_temperature)
    epsilon=zeros(Float64,nlong)
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
        try
            _,M = Wolff!(IS,steps=mlen,save_interval=1)
            M = transpose(M)
            s=reshape(IS.σ,L^2)  
            #matriz[:,i].=s      #agrega los valores a la primera columna de la matriz
            
            r[i],C = space_correlation(binning, reshape(IS.σ,L^2), connected=true,normalized=true)
            #epsilon[i] = correlation_length_r0(r,C)
    catch
            fail += 1
        end
    end
    
    if fail>0 @warn "Failed $(fail)" 
    end
    #writedlm("long_$(L)_$(nlong)_$(mlen).txt",epsilon)
    return r,C
end
