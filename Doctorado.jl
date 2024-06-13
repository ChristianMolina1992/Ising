##1)PROGRAMA PARA CALCULAR LAS MAGNETIZACIONES Y SUSCEPTIBILIDADES, QUE USABA LUEGO PARA VER EL SCALING Y PODER DONDE CAÍA EL MÁXIMO Y USAR ESA TEMPERATURA, DONDE USABA EL ALGORITMO DE METRÓPOLIS PORQUE CON WOLFF NO SE OBSERVAN PÍCOS EN LA SUSCEPTIBILDIAD

using JLD
using LatticeModels
using DelimitedFiles
using Statistics
using BioStatPhys

# Para crear las configuraciones de Wolff

Tc=Ising_SQ_critical_temperature
function config(;T,L,nterm,n)
    IS = Ising(SQLattice_periodic, L, L, ordered=true)

    set_temperature!(IS, T)  #poner el T adecuado, no siempre Tc
    Metropolis!(IS, steps=nterm)
    e,m = Metropolis!(IS,steps=n,save_interval=10)  #calcula la energía y magnetización
    
    #Guardar la configuración en un archivo de texto con formato personalizado, txt y jld
    #ruta = "/home/cmolina/Documentos/Prueba/L=$(L)/Prueba/SQconfig_L$(L)_Temp$(T)_ndata$(n).txt"  
    #file = "/home/cmolina/Documentos/Prueba/L=$(L)/Prueba/SQconfig_L$(L)_Temp$(T)_ndata$(n).jld"
    
    #@save file IS.σ
    #writedlm(ruta, IS.σ)

    mag = Statistics.mean(abs.(m))   
    sus = L * L * (Statistics.var(abs.(m)) / T)
    return mag,sus

end

##2)CONFIGURACIONES YA TERMALIZADAS MEDIANTE WOLFF QUE USARÉ LUEGO PARA HACER METRÓPOLIS
using Statistics

Tc=Ising_SQ_critical_temperature
Temperaturas_50 = collect(1.0:0.01:1.4)*Tc    #el rango que uso cercano a T crítico para ver la susceptibilidad y magnetización vs T y poder ver que está termalizado (cumple el scaling)

function sus_mag(;L)
    magnetizacion = Float64[]
    susceptibility= Float64[]
    
    for T in Temperaturas_50
        mag,sus = config(T=T,L=L,nterm=5000,n=135000)
        push!(magnetizacion,mag)
        push!(susceptibility,sus)
    end
    return magnetizacion,susceptibility
end

mag_50, sus_50 = sus_mag(L=50)


##3)A LA CONFIGURACIÓN OBTENIDA DEL PROGRAMA ANTERIOR PARA LA TEMPERATURA CORRESPONDIENTE Y YA TERMALIZADA, CORRÍA 100000 PASOS MÁS CON LA CONFIGURACIÓN DE WOLFF

using JLD
using DelimitedFiles
using LatticeModels
using BioStatPhys

T = 1.1447916444863286*Tc

function config_wolff(;L,ndata,n)
    IS = Ising(SQLattice_periodic, L, L, ordered=true)
    conf = load("/home/cmolina/Documentos/Prueba/L=$(L)/Prueba/SQconfig_L$(L)_Tc_ndata$(ndata).jld")
    IS.σ = conf["IS.σ"]
    set_temperature!(IS, T)
    set_energy_mag!(IS)
    
    #Rutas donde voy a guardar en formato jld y txt
    ruta = "/home/cmolina/Documentos/Prueba/L=$(L)/Prueba/SQconfig_L$(L)_ndata$(n+ndata).txt"
    file="/home/cmolina/Documentos/Prueba/L=$(L)/Prueba/SQconfig_L$(L)_ndata$(n+ndata).jld"
    
    Wolff!(IS, steps=n, save_interval=1)
    
    #Para guardar 
    
    @save file IS.σ
    writedlm(ruta, IS.σ)    
end

##4) AHORA LO QUE HACÍA ERA, A PARTIR DE LA CONFIGURACIÓN DE WOLFF OBTENIDA EN EL PROGRAMA ANTERIOR, EVOLUCIONABA CON EL ALGORITMO DE WOLFF Y GUARDABA LA CONFIGURACIÓN CADA 1000 PASOS, PARA QUE SEAN CONFIGURACIONES TOTALMENTE DESCORRELACIONADAS. 

using JLD
using LatticeModels
using DelimitedFiles
using BioStatPhys

Tc = Ising_SQ_critical_temperature


function distintas_config_wolff(;L,n,ndata)
    IS = Ising(SQLattice_periodic, L, L, ordered=true)
    conf = load("/home/cmolina/Documentos/Prueba/L=100/Prueba/SQconfig_L$(L)_ndata$(ndata).jld")
    IS.σ = conf["IS.σ"]
    set_temperature!(IS, T)  #poner la correspondiente temperatura
    set_energy_mag!(IS)

    # Ruta base donde se guardarán los archivos
    ruta_base = "/home/cmolina/Documentos/Prueba/L=100/Prueba/SQconfig_L$(L)_pasos"

    # Para guardar las distintas configuraciones cada 1000 pasos
    for i in 0:1000:n
        Wolff!(IS, steps=1000, save_interval=1)
        pasos = ndata + i  # Calcula el número de pasos actual
        ruta = string(ruta_base, pasos, ".txt")  # Construye el nombre del archivo
        writedlm(ruta, IS.σ)
    end
end

##5) EL PASO 5 SERÍA CORRER EL PROGRAMA DE LETI, CON EL CUAL A PARTIR DEL ALGORITMO DE METROPOLIS, IBA OBTENIENDO PARA CADA  CONFIGURACIÓN DE WOLFF (ES DECIR, MI R_0) UNA DISTRIBUCION DE TAUS.

##6) ESTE PROGRAMA ES PARA CALCULAR LAS CORRELACIONES TEMPORALES Y SU RESPECTIVOS TIEMPO DE RELAJACIÓN

using Printf
using JLD
using DelimitedFiles
using LatticeModels
using BioStatPhys

ruta_base = "/home/cmolina/Documentos/Prueba/L=20/Prueba/Muchas metropolis/"

num_archivos = 30  #Cambiar seǵun la cantidad de metrópolis que haya corrido

num_semillas = 25 # Cambiar según la cantidad de semillas que tengas

# Crear un vector de vectores para almacenar las correlaciones de todas las semillas
correlaciones_totales = Vector{Vector{Float64}}[]

# Crear un vector de vectores para almacenar los tiempos de relajación de todas las semillas
tiempos_totales_20 = Vector{Float64}[]

for semilla in 1:num_semillas
    corr = Vector{Float64}[]  # Reiniciar el vector de correlaciones para cada semilla
    tiempos_semilla = Float64[]  # Reiniciar el vector de tiempos de relajación para cada semilla 
    
    for archivo in 1:num_archivos
        # Construir la ruta completa del archivo
        ruta_completa = joinpath(ruta_base, "SQ_L0020_seed$semilla", "Mag-SQconf_L0020_seed$semilla"* "_$(@sprintf("%04d", archivo))")  #acá construyo el nombre completo
        data = readdlm(ruta_completa, header=false, skipstart=4)  #solo me interesa la a partir de la fila 4, por eso el skipstart=4

        M = data[:,2] #solo me interesa la magnetización, por eso elijo esta columna 2
        C = time_correlation(M, connected=true, normalized=true) # Promediando sobre t0
        tiempo_relajacion = correlation_time_spectral(C, 1)  # Calcular el tiempo de relajación
        push!(tiempos_semilla, tiempo_relajacion)
        push!(corr, C)
    end
    push!(tiempos_totales_20, tiempos_semilla)  # Agregar los tiempos de relajación de la semilla al vector total
    push!(correlaciones_totales, corr)   # Agregar la correlación total de la semilla al vector total
end

# Alternativamente, si prefieres concatenar los vectores en lugar de transponer las matrices:
correlaciones_concatenadas_20 = vcat(correlaciones_totales...)

# Escribir los datos concatenados en el archivo de texto
archivo_salida_concatenado = "/home/cmolina/Documentos/Prueba/L=20/Prueba/Muchas metropolis/correlaciones_20_concatenado.txt"

writedlm(archivo_salida_concatenado, correlaciones_concatenadas_20, ' ')
writedlm("/home/cmolina/Documentos/Prueba/L=20/Prueba/Muchas metropolis/tiempos_20.txt",tiempos_totales_20)  #guardo los tiempos de relajación


##7) Codigo para calcular el promedio de las correlaciones, sin usar el readdlm y poder así preservar la memoria

using Statistics

# Definir el nombre del archivo
archivo = "/home/cmolina/Documentos/Prueba/L=100/Prueba/correlaciones_100_concatenado.txt"

# Definir una función para calcular el promedio por columna
function calcular_promedio_columnas(archivo)
    # Crear un vector para almacenar la suma de cada columna
    suma_columnas = zeros(Float64, 1500000)   #acá tengo que tener la longitud de los vectores, por eso con el programa de abajo calculo la longitud del primer vector
    
    # Contador para mantener el número de líneas
    cantidad_filas = 0
    
    # Abrir el archivo
    open(archivo) do file
        # Leer cada línea del archivo
        for linea in eachline(file)
            # Incrementar el contador de líneas
            cantidad_filas += 1
            
            # Dividir la línea en elementos
            elementos = parse.(Float64, split(linea))[1:1500000]
            
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

##7.a)Para calcular cuantos elementos tienen los vectores

# Definir el nombre del archivo
archivo = "/home/cmolina/Documentos/Prueba/L=100/Prueba/correlaciones_100_concatenado.txt"

# Abrir el archivo
open(archivo) do file
    # Leer la primera línea del archivo
    primera_linea = readline(file)
    
    # Contar el número de elementos en la primera línea
    longitud_primer_vector = count(x -> x == ' ', primera_linea) + 1
    
    # Imprimir la longitud del primer vector
    println("La longitud del primer vector es: $longitud_primer_vector")
end

##7.b)Calculo el promedio de las correlaciones y usando el programa 7 y despues lo guardo para realizar el gráfico de las correlaciones medias temporales

using DelimitedFiles

promedio_columnas_100 = calcular_promedio_columnas(archivo)

writedlm("/home/cmolina/Documentos/Prueba/L=100/Prueba/promedio.txt",promedio_columnas_100)

##8)Calculo las distribuciones de los tiempos
using BioStatPhys

@time function histograma_new(;L, bins)
    time = readdlm("/home/cmolina/Documentos/Prueba/L=$(L)/Prueba/tiempos_$(L).txt")
    min,max=extrema(time)
    his = Histogram(bins, max = max, min = min)
        for tiempo in time
            push!(his,tiempo)
        end
    return prob(his)
end
