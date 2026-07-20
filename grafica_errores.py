import pandas as pd
import matplotlib.pyplot as plt

# Leer el archivo CSV generado por C++
df = pd.read_csv('error_vs_tiempo.csv')

plt.figure(figsize=(9, 5), dpi=300)

plt.plot(df['tiempo'], df['error_rms'], label='Error Global RMS ($L_2$)', color='#1d3557', linewidth=2.5)
plt.plot(df['tiempo'], df['error_max'], label='Error Máximo Absoluto ($L_\infty$)', color='#e63946', linestyle='--', linewidth=2.5)

plt.title('Evolución Temporal del Error Numérico (RK4 vs. Solución Analítica)', fontsize=12, fontweight='bold')
plt.xlabel('Tiempo $t$ (s)', fontsize=11)
plt.ylabel('Error de Temperatura (°C)', fontsize=11)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(fontsize=10)
plt.tight_layout()

# Guardar imagen para la diapositiva
plt.savefig('grafica_error_transitorio.png')
print("Gráfica guardada como 'grafica_error_transitorio.png'")
