##1) Como tenìa problema con los transitorios, entonces con este programa agarraba una configuracion y la dejaba evolucionar con wolff, demasiado, como n =70000 o 90000

using JLD
using DelimitedFiles
using LatticeModels
using BioStatPhys

Tc = Ising_SQ_critical_temperature

function config_wolff(;L,ndata,n)
    IS = Ising(SQLattice_periodic, L, L, ordered=true)
    conf = load("/media/cmolina/Datos/Iflysib/Julia/L=20/SQconfig_L$(L)_n$(ndata).jld")
    IS.σ = conf["IS.σ"]
    set_temperature!(IS, Tc)
    set_energy_mag!(IS)
    
    #Rutas donde voy a guardar en formato jld y txt
    ruta = "/media/cmolina/Datos/Iflysib/Julia/L=$(L)/SQconfig_L$(L)_n$(n+ndata).txt"
    file="/media/cmolina/Datos/Iflysib/Julia/L=$(L)/SQconfig_L$(L)_n$(n+ndata).jld"
    
    Wolff!(IS, steps=n, save_interval=1)
    
    #Para guardar 
    
    @save file IS.σ
    writedlm(ruta, IS.σ)    
end




##2) A la configuracion que obtengo, lo que hacía era evolucionar el sistema y cada mil pasos obtenía una nueva configuracion, es decir un nuevo r_0 que serían el total de semillas que tendría



using JLD
using LatticeModels
using DelimitedFiles
using BioStatPhys

Tc = Ising_SQ_critical_temperature


function distintas_config_wolff(;L,n,ndata)
    IS = Ising(SQLattice_periodic, L, L, ordered=true)
    conf = load("/media/cmolina/Datos/Iflysib/Julia/L=20/SQconfig_L$(L)_n$(ndata).jld")
    IS.σ = conf["IS.σ"]
    set_temperature!(IS, Tc)
    set_energy_mag!(IS)

    # Ruta base donde se guardarán los archivos
    ruta_base = "/media/cmolina/Datos/Iflysib/Julia/L=20/SQconfig_L$(L)_pasos"

    # Para guardar las distintas configuraciones cada 1000 pasos
    for i in 0:1000:n
        Wolff!(IS, steps=1000, save_interval=1)
        pasos = ndata + i  # Calcula el número de pasos actual
        ruta = string(ruta_base, pasos, ".txt")  # Construye el nombre del arch>
        writedlm(ruta, IS.σ)
    end
end




#3) El paso tres sería, ya habiendo obtenido todas las semillas anteriores, correr el programa de Leti para obtener por cada semilla, la cantidad de metròpolis distintas que me gustarìa obtener



#4) Acá podía obtener la correlación de cada semilla, es decir, leía la magnetización de cada semilla y cada metròpolis y obtenía un vector donde cada fila se correspondía a una correlación distinta.

using Printf
using JLD
using DelimitedFiles
using LatticeModels
using BioStatPhys

ruta_base = "/media/cmolina/Datos/Iflysib/Julia/L=200/"
#ruta_base = "/media/cmolina/Datos/"

num_archivos = 50

num_semillas = 3 # Cambiar según la cantidad de semillas que tengas

# Crear un vector de vectores para almacenar las correlaciones de todas las semillas
correlaciones_totales = Vector{Vector{Float64}}[]

for semilla in 1:num_semillas
    corr = Vector{Float64}[]  # Reiniciar el vector de correlaciones para cada semilla
    
    for archivo in 1:num_archivos
        # Construir la ruta completa del archivo
        ruta_completa = joinpath(ruta_base, "SQ_L0200_seed$semilla", "Mag-SQconf_L0200_seed$semilla"* "_$(@sprintf("%04d", archivo))")
        data = readdlm(ruta_completa, header=false, skipstart=4)

        M = data[50000:end,2]
        C = time_correlation(M, connected=true, normalized=true) # Promediando sobre t0
        push!(corr, C)
    end
    
    push!(correlaciones_totales, corr)
end

# Alternativamente, si prefieres concatenar los vectores en lugar de transponer las matrices:
correlaciones_concatenadas_200 = vcat(correlaciones_totales...)

# Escribir los datos concatenados en el archivo de texto
archivo_salida_concatenado = "/media/cmolina/Datos/Iflysib/Julia/L=200/correlaciones_200_concatenado.txt"
#archivo_salida_concatenado = "/media/cmolina/Datos/correlaciones_20_concatenado.txt"

writedlm(archivo_salida_concatenado, correlaciones_concatenadas_200, ' ')



#5) Con este programa hacía casi lo miso que antes pero en realidad lo usaba solo para calcular los taus, es decir,, los tiempo de relajaciòn promediando sobre t0

using Printf
using JLD
using DelimitedFiles
using LatticeModels
using BioStatPhys

ruta_base = "/media/cmolina/Datos/Iflysib/Julia/L=200/"
#ruta_base = "/media/cmolina/Datos/"
num_archivos = 50
num_semillas = 3 # Cambiar según la cantidad de semillas que tengas

# Crear un vector de vectores para almacenar los tiempos de relajación de todas las semillas
tiempos_totales_200 = Vector{Float64}[]

for semilla in 1:num_semillas
    tiempos_semilla = Float64[]  # Reiniciar el vector de tiempos de relajación para cada semilla
    
    for archivo in 1:num_archivos
        # Construir la ruta completa del archivo
        ruta_completa = joinpath(ruta_base, "SQ_L0200_seed$semilla", "Mag-SQconf_L0200_seed$semilla"* "_$(@sprintf("%04d", archivo))")
        data = readdlm(ruta_completa, header=false, skipstart=4)

        M = data[50000:end,2]
        C = time_correlation(M, connected=true, normalized=true) # Calcular la correlación
        tiempo_relajacion = correlation_time_spectral(C, 1)  # Calcular el tiempo de relajación
        push!(tiempos_semilla, tiempo_relajacion)
    end
    
    push!(tiempos_totales_200, tiempos_semilla)  # Agregar los tiempos de relajación de la semilla al vector total
end 


