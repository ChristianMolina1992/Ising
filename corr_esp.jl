##Calculo la longitud de correlacion 

using JLD
using LatticeModels
using BioStatPhys
using DelimitedFiles

function compute_long(;L,nlong,mlen,T)
    IS = Ising(SQLattice_periodic,L,L)
    conf = load("/Users/christian/SQconf_Tc_L$(L).jld","IS.σ")
    IS.σ .= conf
    set_energy_mag!(IS)
    set_temperature!(IS,Ising_SQ_critical_temperature)
    epsilon=zeros(Float64,nlong)
    fail = 0
    
    i = 1
    pos = zeros(L, 2)
    for x = 1:L^0.5, y = 1:L^0.5
        pos[i, 1] = x
        pos[i, 2] = y
        i += 1
    end
     binning = distance_binning(pos, 1.0, rmin=0.0)
    
    for i ∈ eachindex(epsilon)
        try
            _,M = Metropolis!(IS,steps=mlen,save_interval=1)
            M = transpose(M)
            #C=BioStatPhys.time_correlation_tw_direct(M,connected=true,i0=1,Xmean=zeros(size(M)),normalized=true)
            # C=time_correlation(M,connected=true,normalized=true,i0=1)
            r,C=space_correlation(binning, load("/Users/christian/SQconf_Tc_L$(L).jld","IS.σ"), connected=true,normalized=true)
            epsilon[i] = correlation_length_r0(r,C)
            #times[i]=correlation_time_spectral(C,1)
        catch
            fail += 1
        end
    end
    if fail>0 @warn "Failed $(fail)" end

    writedlm("long_$(L)_$(nlong)_$(mlen).txt",epsilon)
    return epsilon
end

#Con esto corro el programa @time compute_long(L=20,nlong=1000,mlen=500,T=Ising_SQ_critical_temperature), pero me dan toda las longtudes iguales


##Para las distribuciones y poder graficarlo
using BioStatPhys
@time function histograma_new_long(L,nlong,mlen,bins)
long = readdlm("long_$(L)_$(mlong)_$(mlen).txt",  '\t', '\n')
    #remove(vec(readdlm("times_$(20)_$(500)_$(4420).txt")),0)
min,max=extrema(long)
    his = Histogram(bins, max = max, min = min)
    for longitud in long
        push!(his,longitud)
    end
    return prob(his)
end
