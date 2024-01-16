using JLD
using LatticeModels
using DelimitedFiles

# Para crear las configuraciones de Wolff

L = 50

n = 5000
nterm = 5000
ndata = n + nterm
Magnetizacion = []
# Tc=Ising_SQ_critical_temperature

function configuracion(;T,L,nterm,n)
    for i in 1:ndata
        IS = Ising(SQLattice_periodic, L, L, ordered=true)
        set_temperature!(IS, T)
        Wolff!(IS, steps=ndata)
        
        # Guardar la configuración en un archivo de texto con formato personalizado
        ruta = "/home/cmolina/Documentos/Julia/SQconfig_L$(L)_Temp$(T)_ndata$(i).txt"
        writedlm(ruta, IS.σ)
        mag=IS.M
        push!(Magnetizacion,mag)
    end
end

#Acordate de cmo calcular la susceptibilidad y poner la temperatura critica, y ademas normalizar la temperatura.
