using Printf
using JLD
using DelimitedFiles
using LatticeModels
using BioStatPhys

function procesar_datos_en_bloques(; L, num_archivos, num_semillas, archivos_por_bloque=500)
    ruta_base = "/home/cmolina/L=$L/"

    # Contador de errores en los tiempos de relajación
    contador_errores = 0

    for bloque in 0:Int(num_semillas ÷ archivos_por_bloque)
        # Crear un archivo temporal para cada bloque de 500 semillas
        ruta_salida = joinpath(ruta_base, "tiempos_$(L)_bloque_$(bloque+1).txt")
        tiempos_bloque = Float64[]  # Reiniciar los tiempos para cada bloque

        # Procesar hasta archivos_por_bloque semillas, o menos si es el último bloque
        for semilla in (bloque * archivos_por_bloque + 1):min((bloque + 1) * archivos_por_bloque, num_semillas)
            tiempos_semilla = Float64[]  # Reiniciar el vector de tiempos para cada semilla

            for archivo in 1:num_archivos
                ruta_completa = joinpath(ruta_base, "SQ_L0050_seed$semilla", "Mag-SQconf_L0050_seed$semilla"* "_$(@sprintf("%04d", archivo))")
                
                # Verificar si el archivo existe y no está vacío
                if !isfile(ruta_completa) || filesize(ruta_completa) == 0
                    @warn "Archivo no encontrado o vacío: $ruta_completa"
                    contador_errores += 1
                    continue
                end

                try
                    data = readdlm(ruta_completa, header=false, skipstart=4)

                    if size(data, 1) == 0
                        @warn "Archivo vacío tras leer datos: $ruta_completa"
                        contador_errores += 1
                        continue
                    end

                    M = data[:, 2]  # Solo me interesa la columna de magnetización (columna 2)
                    M = transpose(M)
                    C = BioStatPhys.time_correlation_tw_direct(M, connected=true, i0=1, Xmean=zeros(size(M)), normalized=true)  # Sin promediar sobre t0
                    
                    # Intentar calcular el tiempo de relajación
                    tiempo_relajacion = correlation_time_spectral(C, 1)
                    push!(tiempos_semilla, tiempo_relajacion)

                catch e
                    @warn "Error al procesar archivo $archivo de semilla $semilla: $e"
                    contador_errores += 1
                    continue
                end
            end

            # Agregar los tiempos de la semilla al bloque
            append!(tiempos_bloque, tiempos_semilla)
        end

        # Guardar los tiempos del bloque en un archivo
        writedlm(ruta_salida, tiempos_bloque)
        println("Guardado bloque $(bloque+1) en $ruta_salida")
    end

    # Concatenar todos los bloques al final
    concatenar_archivos_bloques(L, ruta_base, num_semillas, archivos_por_bloque)

    # Imprimir el número total de errores
    println("Número total de errores en el cálculo de tiempos de relajación: $contador_errores")
end

# Función para concatenar todos los archivos en un solo archivo
function concatenar_archivos_bloques(L, ruta_base, num_semillas, archivos_por_bloque)
    archivo_salida_final = joinpath(ruta_base, "tiempos_$(L)_completo.txt")
    open(archivo_salida_final, "w") do io
        for bloque in 0:Int(num_semillas ÷ archivos_por_bloque)
            ruta_bloque = joinpath(ruta_base, "tiempos_$(L)_bloque_$(bloque+1).txt")
            if isfile(ruta_bloque)
                data_bloque = readdlm(ruta_bloque)
                writedlm(io, data_bloque)
            end
        end
    end
    println("Todos los bloques concatenados en $archivo_salida_final")
end

# Uso del programa
L = 50
num_archivos = 1  # Sigues trabajando con un solo archivo por semilla
num_semillas = 6839  # El total de semillas que quieres procesar

# Procesar en bloques de 500 semillas
procesar_datos_en_bloques(L=L, num_archivos=num_archivos, num_semillas=num_semillas)


