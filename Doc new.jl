#Para calcular las correlaciones para todas las semillas 

using Printf
using JLD
using DelimitedFiles
using LatticeModels
using BioStatPhys

ruta_base = "/media/cmolina/Datos/Iflysib/Julia/L=200/"
num_archivos = 27
num_semillas = 1 # Cambiar según la cantidad de semillas que tengas

# Crear un vector de vectores para almacenar las correlaciones de todas las semillas
correlaciones_totales = Vector{Vector{Float64}}[]

for semilla in 1:num_semillas
    corr = Vector{Float64}[]  # Reiniciar el vector de correlaciones para cada semilla
    
    for archivo in 1:num_archivos
        # Construir la ruta completa del archivo
        ruta_completa = joinpath(ruta_base, "SQ_L0200_seed$semilla", "Mag-SQconf_L0200_seed$semilla"* "_$(@sprintf("%04d", archivo))")
        data = readdlm(ruta_completa, header=false, skipstart=4)

        M = data[:, 2]
        C = time_correlation(M, connected=true, normalized=true) # Promediando sobre t0
        push!(corr, C)
    end
    
    push!(correlaciones_totales, corr)
end

# Alternativamente, si prefieres concatenar los vectores en lugar de transponer las matrices:
correlaciones_concatenadas_200 = vcat(correlaciones_totales...)

# Escribir los datos concatenados en el archivo de texto
archivo_salida_concatenado = "/media/cmolina/Datos/Iflysib/correlaciones_200_concatenado.txt"
writedlm(archivo_salida_concatenado, correlaciones_concatenadas_200, ' ')






#Este es para calcular los taus 

using Printf
using JLD
using DelimitedFiles
using LatticeModels
using BioStatPhys

ruta_base = "/media/cmolina/Datos/Iflysib/Julia/L=200/"
num_archivos = 27
num_semillas = 1 # Cambiar según la cantidad de semillas que tengas

# Crear un vector de vectores para almacenar los tiempos de relajación de todas las semillas
tiempos_totales_200 = Vector{Float64}[]

for semilla in 1:num_semillas
    tiempos_semilla = Float64[]  # Reiniciar el vector de tiempos de relajación para cada semilla
    
    for archivo in 1:num_archivos
        # Construir la ruta completa del archivo
        ruta_completa = joinpath(ruta_base, "SQ_L0200_seed$semilla", "Mag-SQconf_L0200_seed$semilla"* "_$(@sprintf("%04d", archivo))")
        data = readdlm(ruta_completa, header=false, skipstart=4)

        M = data[:, 2]
        C = time_correlation(M, connected=true, normalized=true) # Calcular la correlación
        tiempo_relajacion = correlation_time_spectral(C, 1)  # Calcular el tiempo de relajación
        push!(tiempos_semilla, tiempo_relajacion)
    end
    
    push!(tiempos_totales_200, tiempos_semilla)  # Agregar los tiempos de relajación de la semilla al vector total
end 





##Codigo para calcular el promedio de las correlaciones, sin usar el readdlm y preservar la memoria

using Statistics

# Definir el nombre del archivo
archivo = "/media/cmolina/Datos/Iflysib/correlaciones_200_concatenado.txt"

# Definir una función para calcular el promedio por columna
function calcular_promedio_columnas(archivo)
    # Crear un vector para almacenar la suma de cada columna
    suma_columnas = zeros(Float64, 750000)
    
    # Contador para mantener el número de líneas
    cantidad_filas = 0
    
    # Abrir el archivo
    open(archivo) do file
        # Leer cada línea del archivo
        for linea in eachline(file)
            # Incrementar el contador de líneas
            cantidad_filas += 1
            
            # Dividir la línea en elementos
            elementos = parse.(Float64, split(linea))
            
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

# Calcular los promedios de las columnas
promedio_columnas_200 = calcular_promedio_columnas(archivo)




##Para calcular cuantos elementos tienen los vectores

# Definir el nombre del archivo
archivo = "/media/cmolina/Datos/Iflysib/correlaciones_200_concatenado.txt"

# Abrir el archivo
open(archivo) do file
    # Leer la primera línea del archivo
    primera_linea = readline(file)
    
    # Contar el número de elementos en la primera línea
    longitud_primer_vector = count(x -> x == ' ', primera_linea) + 1
    
    # Imprimir la longitud del primer vector
    println("La longitud del primer vector es: $longitud_primer_vector")
end
