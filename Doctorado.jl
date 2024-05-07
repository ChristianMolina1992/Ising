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


#6)##Codigo para calcular el promedio de las correlaciones, sin usar el readdlm y preservar la memoria

using Statistics

# Definir el nombre del archivo
archivo = "/media/cmolina/Datos/Iflysib/Julia/L=20/correlaciones_20_concatenado.txt"

# Definir una función para calcular el promedio por columna
function calcular_promedio_columnas(archivo)
    # Crear un vector para almacenar la suma de cada columna
    suma_columnas = zeros(Float64, 499000)
    
    # Contador para mantener el número de líneas
    cantidad_filas = 0
    
    # Abrir el archivo
    open(archivo) do file
        # Leer cada línea del archivo
        for linea in eachline(file)
            # Incrementar el contador de líneas
            cantidad_filas += 1
            
            # Dividir la línea en elementos
            elementos = parse.(Float64, split(linea))[1:499000]
            
            # Sumar cada elemento a la suma correspondiente de la columna
            for (indice, valor) in enumerate(elementos)
                suma_columnas[indice] += valor
            end
        end
    end
    
    # Calcular el promedio de cada columna
    promedio_columnas = suma_columnas / cantidad_filas
    
    return promedio_columnas
end



#Y con esto guardaba los taus
using DelimitedFiles

# Ruta para guardar el archivo
ruta_directorio = "/media/cmolina/Datos/"
ruta_guardado = joinpath(ruta_directorio, "promedio_columnas_200.txt")

# Crear el directorio si no existe
if !isdir(ruta_directorio)
    mkdir(ruta_directorio)
end

# Guardar los datos en el archivo de texto
writedlm(ruta_guardado, promedio_columnas_200, ' ')


# Con esto calculaba la cantidad de vectores que tenía la primera fila para usarlo arriba
promedio_columnas_20 = calcular_promedio_columnas(archivo)





##Para calcular cuantos elementos tienen los vectores

# Definir el nombre del archivo
archivo = "/media/cmolina/Datos/Iflysib/Julia/L=20/correlaciones_20_concatenado.txt"

# Abrir el archivo
open(archivo) do file
    # Leer la primera línea del archivo
    primera_linea = readline(file)
    
    # Contar el número de elementos en la primera línea
    longitud_primer_vector = count(x -> x == ' ', primera_linea) + 1
    
    # Imprimir la longitud del primer vector
    println("La longitud del primer vector es: $longitud_primer_vector")
end



## Programa para cambiar nombres

# Función para cambiar el nombre de los archivos en un directorio
function cambiar_nombres(directorio::String, prefijo_original::String, prefijo_nuevo::String)
    # Obtener lista de archivos en el directorio
    archivos = readdir("/media/cmolina/Datos/Iflysib/SQ_L0100_seed4/")
    
    # Iterar sobre cada archivo
    for archivo in archivos
        # Verificar si el archivo comienza con el prefijo original
        if startswith(archivo, prefijo_original)
            # Construir el nuevo nombre de archivo reemplazando el prefijo
            nuevo_nombre = replace(archivo, prefijo_original => prefijo_nuevo)
            # Renombrar el archivo
            mv(joinpath(directorio, archivo), joinpath(directorio, nuevo_nombre))
        end
    end
end

# Directorio donde se encuentran los archivos
directorio = "/media/cmolina/Datos/Iflysib/SQ_L0100_seed4/"

# Prefijo original y nuevo
prefijo_original = "Mag-SQconf_L0100_seed4"
prefijo_nuevo = "Mag-SQconf_L0100_seed16"

# Llamar a la función para cambiar los nombres de los archivos
cambiar_nombres(directorio, prefijo_original, prefijo_nuevo)

