##CREO LAS DISTINTAS CONFIGURTACIONES CON EL ALGORITMO DE WOLFF

using JLD

Ls=[20,50,100,200,500]
Tc=Ising_SQ_critical_temperature
for L ∈ Ls
    IS=Ising(SQLattice_periodic,L,L,ordered=true)
    set_temperature!(IS,Tc)
    Wolff!(IS,steps=5000)
    file = IS.σ
    writedlm("wolf_$(L).txt",file)
     
    println("Saved L=$L")
end

using BioStatPhys
n_sweep = 1000000
@time function correlacion_temporal_wolf(T,L)
    spin = readdlm("wolf_$(L).txt", Int8)
    m1=Float64[]
    for i=1:n_sweep
        ising2d_ifelse!(spin,1/T,1)
        push!(m1,magnetization_ising2d(spin))                                  
    end
     corr=time_correlation(m1,connected=true,normalized=true)
    time=collect(1:size(corr,1))
    return corr
    #return plot(time, corr, xlabel = "time", ylabel = "C(t)")  
end
