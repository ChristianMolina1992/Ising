#Para graficar y ajustar a una recta

import numpy as np
import matplotlib.pyplot as plt

# Datos
x = np.array([50, 100,200])
y = np.array([444, 9845, 13111])

# Ajuste lineal
m, b = np.polyfit(x, y, 1)

# Calcular el valor de la pendiente
pendiente = m

# Generar puntos para la línea ajustada
x_fit = np.linspace(min(x), max(x), 100)
y_fit = m * x_fit + b

# Graficar los datos y la línea ajustada
plt.figure()
plt.scatter(x, y, label='Datos')
plt.plot(x_fit, y_fit, 'r-', label=f'Ajuste lineal (Pendiente={pendiente:.2f})')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Ajuste lineal')
plt.legend()
plt.grid(True)
plt.show()

# Imprimir el valor de la pendiente
print("La pendiente del ajuste lineal es:", pendiente)
print("orrd",b)




import numpy as np
import matplotlib.pyplot as plt

# Datos
x = np.array([20,50,100])
y = np.array([13688,34004,306310])

# Ajuste lineal
m, b = np.polyfit(np.log(x), np.log(y), 1)

# Calcular el valor de la pendiente
pendiente = m

# Generar puntos para la línea ajustada
x_fit = np.linspace(min(x), max(x), 100)
y_fit = np.exp(m * np.log(x_fit) + b)

# Graficar los datos y la línea ajustada en escala logarítmica
plt.figure()
plt.scatter(x, y, label='Datos')
plt.plot(x_fit, y_fit, 'r-', label='Ajuste lineal')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('L')
plt.ylabel('Tau')
plt.title('Tau vs L en escala logarítmica')

# Establecer marcas en los ejes x e y con los valores de los puntos de datos en escala logarítmica
plt.xticks(x, x)
plt.yticks(y, y)

plt.legend()
plt.grid(True)
plt.show()

# Imprimir el valor de la pendiente
print("La pendiente del ajuste lineal es:", pendiente)


import numpy as np
import matplotlib.pyplot as plt

# cargo los datos
data = np.loadtxt("/home/cmolina/Documentos/Nuevo/L=50/SQ_L0050_seed1/Mag-SQconf_L0050_seed1_0010", skiprows=4)
mag = data[:, 1]

# Defino rango
x_range = np.arange(1, 2500001)[:1000000]
mag_subset = mag[:1000000]

plt.plot(x_range, mag_subset)
plt.xlabel('Tiempo')
plt.ylabel('Magnetización')
plt.title('Magnetización vs tiempo')
plt.show()

