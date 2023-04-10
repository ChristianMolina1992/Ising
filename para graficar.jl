#Corro esto p_20 = @time compute_times(L=20, ntimes=2000, mlen=40000)

#Leo esto using DelimitedFiles 
p_50=readdlm("times_$(50)_$(2000)_$(250000).txt")


#Remuevo los ceros que fallan
using Statistics 
function remove(g,a)
    deleteat!(g, findall(x->x==a, g))
    return g
end

#Lo corro al programa de arriba
ep_100 =remove(vec(p_100),0)


#Calculo las funciones de distribucion
using BioStatPhys
@time function histograma_new(L,ntimes,mlen,bins)
time = readdlm("times_$(L)_$(ntimes)_$(mlen).txt",  '\t', '\n')
    #remove(vec(readdlm("times_$(20)_$(500)_$(4420).txt")),0)
min,max=extrema(time)
    his = Histogram(bins, max = max, min = min)
    for tiempo in time
        push!(his,tiempo)
    end
    return prob(his)
end

# Con esto hago todo lo de arriba para obtener directamente lo que queiro graficar
function arre(L,ntimes,mlen,bins)
    p_L=readdlm("times_$(L)_$(ntimes)_$(mlen).txt")
    ep_L =remove(vec(p_L),0)
    t_L,ep_L = histograma_new(L,ntimes,mlen,bins)
    return x_L,y_L = log.(t_L),log.(ep_L)
end

#Con esto corro lo que quiero graficar y obtengo los ejes
x_50,y_50=arre(50,2000,250000,20)


#Grafico
using Plots
SR1 = plot(x_20,y_20,label="L=20",marker=:circle, line=:solid)
SR2 = plot(SR1,x_50,y_50,label="L=50",marker=:circle, line=:solid)
plot!(SR2,x_100,y_100,label="L=100",ylabel = "Densidad de probabilidad en log",xlabel = "Tiempos de relajacion en log",marker=:circle, line=:solid)



#Para graficar las distribuciones como el paper, hago lo siguiente

#Esto seria el eje Y
function dis(L,ntimes,mlen,bins)
    X_L,Y_L=histograma_new(L,ntimes,mlen,bins)
    return X_L.*Y_L
end

#Esto es el eje X

function ejex(L,ntimes,mlen,bins)
    X_L,Y_L=histograma_new(L,ntimes,mlen,bins)
    return X_L/L^2.21
end

#Con esto lo corro
X_50= ejex(50,2000,250000,20)


#Por último, grafico
G1 = plot(X_20,Y_20,label="L=20",marker=:circle, line=:solid)
G2 = plot(G1,X_50,Y_50,label="L=50",marker=:circle, line=:solid)
plot!(G2,X_100,Y_100,label="L=100",ylabel = "P(t)*t",xlabel = "t/L^2.21",marker=:circle, line=:solid)



