using LatticeModels
IS = Ising(SQLattice_periodic,50,50;T=2.7)
println("Energy $(IS.E), magnetization $(IS.M)")

##CREO LAS CONFIGURACIONES

using JLD

Ls=[20,50,100,200,500]
Tc=Ising_SQ_critical_temperature
for L ∈ Ls
    IS=Ising(SQLattice_periodic,L,L,ordered=true)
    set_temperature!(IS,Tc)
    Wolff!(IS,steps=5000)
    file="/Users/christian/SQconf_Tc_L$(L).jld"
    @save file IS.σ
    println("Saved L=$L")
end

##CALCULO LOS TAU

using JLD
using LatticeModels
using BioStatPhys

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
            C=BioStatPhys.time_correlation_tw_direct(M,connected=true,i0=1,
                                                     Xmean=zeros(size(M)),normalized=true)
            # C=time_correlation(M,connected=true,normalized=true,i0=1)
           
            times[i]=correlation_time_spectral(C,1)
        catch
            fail += 1
        end
    end
    if fail>0 @warn "Failed $(fail)" end
    return times
end

##PARA CONTAR LOS CEROS Y VERIFICAR QUE COINCIDE CON LOS TAU FALLIDOS

B = [(0, count(==(0), p_20))]


##PARA ELIMINAR LOS CEROS Y DARME EL VALOR MEDIO

using Statistics 
function remove!(g,a)
    deleteat!(g, findall(x->x==a, g))
    return Statistics.mean(g)
end
