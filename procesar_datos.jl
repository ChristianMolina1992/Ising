include("corr.jl")

# Define los parámetros aquí
L = 50
num_archivos = 1  # Sigues trabajando con un solo archivo por semilla
num_semillas = 100  # El total de semillas que quieres procesar
#archivo_inicio = 1
#archivo_fin = 1
#append_mode = false  # No necesitamos "append" porque solo guardamos todo de una vez

# Procesar todas las semillas en una sola ejecución
println("Procesando todas las semillas (1 a 1000)...")
#procesar_datos(L=L, num_archivos=num_archivos, num_semillas=num_semillas, archivo_inicio=archivo_inicio, archivo_fin=archivo_fin, append_mode=append_mode)
procesar_datos(L=L,num_archivos=num_archivos,num_semillas=num_semillas)
println("Procesamiento completo.")

