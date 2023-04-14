
##EL PROBLEMA QUE TENGO ES QUE ME DA MUY DISTINTO LOS GRAFICOS CADA VEZ QUE CORRO EL PROGRAMA, LE PONGO MLEN=4000 PARA ASEGURARME QUE TERMALIZE PERO ME SIGUE DANDO DISTINTO


using JLD
using LatticeModels
using BioStatPhys
using DelimitedFiles
using Plots 

function prueba_corr_rvs_C(;L,mlen,longitud)
    IS = Ising(SQLattice_periodic,L,L)
    conf = load("/Users/christian/SQconf_Tc_L$(L).jld","IS.σ")
    IS.σ = conf
    matriz=zeros((L^2),longitud)  #matriz de ceros que iré rellenando, de filas L^2 y columnas longitud, esta última indica las distintas configuraciones que iré agregando
    set_energy_mag!(IS)
    set_temperature!(IS,Ising_SQ_critical_temperature)
    #epsilon=zeros(Float64,nlong)
    fail = 0
    
    i = 1
    pos = zeros(L^2, 2)
    for x = 1:L, y = 1:L
        pos[i, 1] = x
        pos[i, 2] = y
        i += 1
    end
     binning = distance_binning(pos, 1.0, rmin=0.0) #Con esto defino la función binning
    
    for i in 1:longitud     #esta funcion me calcula wolff definiendo la cantidad de pasos con mlen y luego la longitud me indicara cuantas configuraciones creo
        try
            _,M = Wolff!(IS,steps=mlen,save_interval=1)
            M = transpose(M)
            s=reshape(IS.σ,L^2)  
            matriz[:,i].=s      #agrega los valores a la primera columna de la matriz, a medida que aumente la longitud, le iré agregando los valores a las columnas
            
            #C=BioStatPhys.time_correlation_tw_direct(M,connected=true,i0=1,Xmean=zeros(size(M)),normalized=true)
            # C=time_correlation(M,connected=true,normalized=true,i0=1)
            #r,C = space_correlation(binning, reshape(IS.σ,L^2), connected=true,normalized=true)
            #epsilon[i] = correlation_length_r0(r,C)
            #times[i]=correlation_time_spectral(C,1)
        catch
            fail += 1
        end
    end
    
    if fail>0 @warn "Failed $(fail)" 
    end
    X=reshape(matriz,(L^2)*longitud)     
    r,C = space_correlation(binning, X, connected=true,normalized=true)

    #writedlm("long_$(L)_$(nlong)_$(mlen).txt",epsilon)
    return r,C
    #return X
    
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
