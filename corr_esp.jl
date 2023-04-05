###Esto no da, el error esta en el binning, si corro solo el binning no da


using LatticeModels
using BioStatPhys
using DelimitedFiles

# Definir el tamaño de la red y el binning espacial
L = 50
pos = collect([(x,y) for x in 1:L for y in 1:L], L, L)
binning = distance_binning(pos, 1.)

# Cargar los datos de la variable X
X = load("X_$(L).txt")

# Calcular la correlación espacial de X
corr = space_correlation(binning, X, connected=true, normalized=true, Xmean=mean(X))

# Guardar los resultados en un archivo de texto
#corr_esp = "corr_esp_$(L).txt.txt"
#writedlm(aca guardaria las correlaciones)


###Esto parecde que funciona,

using BioStatPhys

# Generar 100 puntos aleatorios en un cubo de lado 10
pos = 10 .* rand(100, 3)

# Calcular la matriz de distancias para todos los pares de puntos
bins = distance_binning(pos, 1.0)

# Obtener los índices de los pares de puntos a una distancia de 2 a 3 unidades
pairs_indices = bins[2.0:3.0]
