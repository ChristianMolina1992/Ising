#Para generar las configuraciones termalizadas

using JLD

Ls=[20,50,100,200,500]
Tc=Ising_SQ_critical_temperature
for L ∈ Ls
    IS=Ising(SQLattice_periodic,L,L,ordered=true)
    set_temperature!(IS,Tc)
    Wolff!(IS,steps=5000)
    file="./configs/SQconf_Tc_L$(L).jld"
    @save file IS.σ
    println("Saved L=$L")
end


#Para calcular los tiempos de relajacion a partir de las configuraciones termalizadas

using JLD
using LatticeModels
using BioStatPhys
using DelimitedFiles

function compute_times(;L,ntimes,mlen,T)
    IS = Ising(SQLattice_periodic,L,L)
    conf = load("/Users/christian/SQconf_Tc_L$(L).jld","IS.σ")
    IS.σ .= conf
    set_energy_mag!(IS)
    set_temperature!(IS,T)
    times=zeros(Float64,ntimes)
    fail = 0
    for i ∈ eachindex(times)
        try
            _,M = Metropolis!(IS,steps=mlen,save_interval=1)
            M = transpose(M)
            C=BioStatPhys.time_correlation_tw_direct(M,connected=true,i0=1,
                                                     Xmean=zeros(size(M)),normalized=true)
            # C=time_correlation(M,connected=true,normalized=true,i0=1)
           
            times[i]=correlation_time_spectral(C,1)
        catch
            fail += 1
        end
    end
    if fail>0 @warn "Failed $(fail)" end

    writedlm("times_$(L)_$(ntimes)_$(mlen)_$(T).txt",times)
    return times
end

#Para calcular las longitudes de correlacion r_0 a partir de las configuraciones termalizadas

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
            
            r,C = space_correlation(binning, reshape(IS.σ,L^2), connected=true,normalized=true)
            epsilon[i] = correlation_length_r0(r,C)
    catch
            fail += 1
        end
    end
    
    if fail>0 @warn "Failed $(fail)" 
    end
    writedlm("long_$(L)_$(nlong)_$(mlen).txt",epsilon)
    return epsilon
end

#Para calcular el histograma de una correlacion, me devuelve el P(r_0)*r_0 (no recuerdo bien el orden)

using BioStatPhys
@time function histograma_new_long(;L,nlong,mlen,bins)
long = readdlm("long_$(L)_$(nlong)_$(mlen).txt",  '\t', '\n')
    min,max=extrema(long)
    his = Histogram(bins, max = max, min = min)
    for longitud in long
        push!(his,longitud)
    end
    return prob(his)
end

#Para calcular los r_0 y los tau de una misma configuracion

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

