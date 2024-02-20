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





using JLD
using LatticeModels
using DelimitedFiles
using Statistics

# Para crear las configuraciones de Wolff


function configuracion(;T,L,nterm,n)
        IS = Ising(SQLattice_periodic, L, L, ordered=true)
        set_temperature!(IS, T)
        Wolff!(IS, steps=nterm)
        e,m = Wolff!(IS,steps=n,save_interval=10)
        # Guardar la configuración en un archivo de texto con formato personalizado
        ruta = "/home/cmolina/Documentos/Julia/Config1.0/SQconfig_L$(L)_Temp$(T)_ndata$(n).txt"
        writedlm(ruta, IS.σ)

        mag=mean(abs.(m))
        sus=L*L*var(m)
        return mag,sus
 end



L = 20
Tc=Ising_SQ_critical_temperature

magnetizacion_20 = Float64[]
susceptibility_20 = Float64[]
for T=0.2*Tc:0.05:1.4*Tc
    mag,sus = configuracion(T=T,L=L,nterm=5000,n=10000)
    push!(magnetizacion_20,mag)
    push!(susceptibility_20,sus)
end

using JLD
using LatticeModels
using Plots

Tc = Ising_SQ_critical_temperature
ejey = [0.0223066,0.02209,0.02315817] .* 100^(7/4)
ejex = (([ 0.9, 1.0,1.1,] .* Tc .- Tc) ./ Tc) * (-100)

plot(ejex, ejey, seriestype=:scatter,label="L100")




##PARA VER SI ESTA BIEN


Tc=Ising_SQ_critical_temperature
L=100
n=10000


function mag()
    IS = Ising(SQLattice_periodic, L, L, ordered=true)
    conf=load("/home/cmolina/Documentos/Julia/SQconfig_L100_Tc_ndata20000.jld")
    IS.σ = conf["IS.σ"]
    set_temperature!(IS, Tc)
    set_energy_mag!(IS)
    E,M = Wolff!(IS,steps=n,save_interval=1)
    return M
end
C = BioStatPhys.time_correlation_tw_direct(transpose(M_100), connected=true, i0=1, Xmean=zeros(size(M)), normalized=true)

correlation_time_spectral(C,1)
