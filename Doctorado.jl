using JLD
using LatticeModels
using DelimitedFiles

# Para crear las configuraciones de Wolff. Con este programa puedo ir guardando las configuraciones como así tambien ir calculando las magnetizaciones por cada configuración e ir guardándolas

L = 50

n = 5000
nterm = 5000
ndata = n + nterm
Magnetizacion = []
Tc=Ising_SQ_critical_temperature

function configuracion(;T,L,nterm,n)
    for i in 1:ndata
        IS = Ising(SQLattice_periodic, L, L, ordered=true)
        set_temperature!(IS, T)
        Wolff!(IS, steps=ndata)
        
        # Guardar la configuración en un archivo de texto con formato personalizado
        ruta = "/home/cmolina/Documentos/Julia/Configuraciones/SQconfig_L$(L)_Temp$(T)_ndata$(i).txt"
        writedlm(ruta, IS.σ)
        #mag=IS.M
        push!(Magnetizacion,M)
        writedlm("/home/cmolina/Documentos/Julia/Configuraciones/magnetizacion$(L).txt",M)
    end
end



##Como me di cuenta que el programa de arriba no me arrojaba las magnetizaciones por spin, entonces realicé otro programa para que me la devolviera y finalmente grafíqué haciendo using Plots 
#plot((1:10000),magnetizaciones,label=false,ylabel="Magnetización",xlabel="Tiempo")

function magnetization_ising2d(s)
    m, n = size(s)
    return sum(s)/(m*n)
end

ruta = "/home/cmolina/Documentos/Julia/"
nombre_en_comun = "SQconfig_L50_Temp2.269185314213022_ndata"

cant_archivos = 1:10000
magnetizaciones = []

for i in cant_archivos
    nombre_archivo = joinpath(ruta, "$nombre_en_comun$i.txt")
    
    try
        # Leer el archivo
        s = readdlm(nombre_archivo)
        
        # Aplicar la función y almacenar el resultado
        magnetizacion = abs(magnetization_ising2d(s))
        
        push!(magnetizaciones, magnetizacion)
        
    catch e
    end
end

