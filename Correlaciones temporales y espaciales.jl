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
    correlacion = []
    fail = 0
    for i ∈ eachindex(times)
        try
            _, M = Metropolis!(IS, steps = mlen, save_interval = 1)
            M = transpose(M)
            C= BioStatPhys.time_correlation_tw_direct(M, connected = true, i0 = 1, Xmean = zeros(size(M)), normalized = true)
            #println("Iteración $i: C = $C")
            times[i]=correlation_time_spectral(C,1)
            push!(correlacion,C)
            #eje=collect(1:size(C,1))
        catch
            fail += 1
        end
    end
    if fail > 0
        @warn "Failed $(fail)"
    end
 return correlacion

end


# Asegurate de que todos los vectores tengan la misma longitud
longitud_referencia = length(y_20[2])
#for vec in y_20
    #if length(vec) != longitud_referencia
        #throw(ArgumentError("Los vectores deben tener la misma longitud"))
    #end
#end

# Inicializa la suma con un vector de ceros
suma = zeros(eltype(y_20[1]), longitud_referencia)

# Calcula la suma de los componentes correspondientes en cada vector
for vec in y_20
    if length(vec) == longitud_referencia
        suma .+= vec
    end
end

# Calcula el promedio dividiendo la suma entre el número de vectores originales
num_vectores = length(y_20)
promedio = suma / (num_vectores-1)



using JLD
using LatticeModels
using BioStatPhys
using DelimitedFiles
using Plots 

function correlacion_espacial(;L,nlong,mlen)
    IS = Ising(SQLattice_periodic,L,L)
    conf = load("/Users/christian/SQconf_Tc_L$(L).jld","IS.σ")
    IS.σ = conf
    set_energy_mag!(IS)
    set_temperature!(IS,Ising_SQ_critical_temperature)
    epsilon=zeros(Float64,nlong)
    correlacion = []
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
            
            r,C = space_correlation(binning, reshape(IS.σ,L^2), connected=true,normalized=true)
            push!(correlacion,C)
            #epsilon[i] = correlation_length_r0(r,C)
    catch
            fail += 1
        end
    end
    
    if fail>0 @warn "Failed $(fail)" 
    end
    
    return correlacion
end

#Para graficar
#@time dato_20_corr_esp = correlacion_espacial(L=20,nlong=2000,mlen=20000)
