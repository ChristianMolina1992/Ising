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

##OBTENGO LA CORRELACION TEMPORAL Y LUEG CON CORRELATION_TIME_SPECTRAL CALCULO EL TIEMPO DE RELAJACION PARA CADA TAMAÑO (AMBOS PROGRAMAS PARA TC)

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

##AHORA CREO LAS CONFIGURACIONES MEDIANTE EL ALGORITMO DE METROPOLIS

using JLD

Ls=[20,50,100,200,500]
Tc=Ising_SQ_critical_temperature
for L ∈ Ls
    IS=Ising(SQLattice_periodic,L,L,ordered=true)
    set_temperature!(IS,Tc)
    Metropolis!(IS,steps=5000)
    file = IS.σ
    writedlm("metropolis_$(L).txt",file)
     
    println("Saved L=$L")
end

##OBTENGO LA CORRELACION TEMPORAL DE LAS CONFIGURACIONES

using BioStatPhys
n_sweep = 1000
@time function correlacion_temporal_metropolis(T,L)
    spin = readdlm("metropolis_$(L).txt", Int8)
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


##POR ÚLTIMO CALCULO LAS FUNCIONES Y LAS DISTRIBUCIONES DE TAU (FALTARÍA GRAFICAR)

n_sweep = 200000
@time function distribuciones_metropolis(L)
    spin = readdlm("metropolis_$(L).txt", Int8) 

    tiempo_relaj = Float64[]
    m1 = Float64[]

    for i=1:n_sweep
        ising2d_ifelse!(spin,1/Ising_SQ_critical_temperature,1)
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
    writedlm("tiempos_metropolis$(L).txt",tiempo_relaj)
    return collect(tiempo_relaj)
end

@time function histograma_metropolis(L)
time = readdlm(readdlm("metropolis_$(L).txt"),  '\t', '\n')
min,max=extrema(time)
    his = Histogram(100, max = max, min = min)
    for tiempo in time
        push!(his,tiempo)
    end
    return prob(his)
end
