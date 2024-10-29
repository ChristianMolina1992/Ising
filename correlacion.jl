using Printf
using JLD
using DelimitedFiles
using LatticeModels
using BioStatPhys

function procesar_datos(; L, num_archivos, num_semillas)
    ruta_base = "/home/cmolina/L=$L/"

    # Crear un vector de vectores para almacenar las correlaciones de todas las semillas
    #correlaciones_totales = Vector{Vector{Float64}}[]

    # Crear un vector de vectores para almacenar los tiempos de relajación de todas las semillas
    tiempos_totales = Vector{Float64}[]

    # Contador de errores en los tiempos de relajación
    contador_errores = 0

    for semilla in 1:num_semillas
       # corr = Vector{Float64}[]  # Reiniciar el vector de correlaciones para cada semilla
        tiempos_semilla = Float64[]  # Reiniciar el vector de tiempos de relajación para cada semilla

        for archivo in 1:num_archivos
            ruta_completa = joinpath(ruta_base, "SQ_L0050_seed$semilla", "Mag-SQconf_L0050_seed$semilla"* "_$(@sprintf("%04d", archivo))")
            data = readdlm(ruta_completa, header=false, skipstart=4)

            M = data[100000:end,2] # Solo me interesa la magnetización, por eso elijo esta columna 2
            M = transpose(M)
            C = BioStatPhys.time_correlation_tw_direct(M, connected=true, i0=1, Xmean=zeros(size(M)), normalized=true) # Sin promediar sobre t0
            
            # Intentar calcular el tiempo de relajación y capturar errores
            try
                tiempo_relajacion = correlation_time_spectral(C, 1)  # Calcular el tiempo de relajación
                push!(tiempos_semilla, tiempo_relajacion)
            catch e
                @warn "Error al calcular tiempo de relajación en semilla $semilla, archivo $archivo: $e"
                contador_errores += 1  # Incrementar el contador de errores
                continue  # Saltar al siguiente archivo si ocurre un error
            end

            #push!(corr, C)
        end

        push!(tiempos_totales, tiempos_semilla)  # Agregar los tiempos de relajación de la semilla al vector total
        #push!(correlaciones_totales, corr)       # Agregar la correlación total de la semilla al vector total
    end

    # Concatenar todas las correlaciones
    #correlaciones_concatenadas = vcat(correlaciones_totales...)

    # Escribir los datos concatenados en el archivo de texto
   # archivo_salida_concatenado = joinpath(ruta_base, "correlaciones_$(L)_concatenado.txt")
    
    # Usar modo append si se especifica
    #if append_mode
     #  open(archivo_salida_concatenado, "a") do io
      #      writedlm(io, correlaciones_concatenadas, ' ')
       # end
    #else
     #   writedlm(archivo_salida_concatenado, correlaciones_concatenadas, ' ')
   # end

    writedlm(joinpath(ruta_base, "tiempos_$(L).txt"), tiempos_totales)  #guardo los tiempos de relajación

    # Imprimir el número total de errores
    println("Número total de errores en el cálculo de tiempos de relajación: $contador_errores")
end
