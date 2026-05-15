import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

os.environ.pop('mplbackend', None)
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", palette="muted")

def cargar_datos_expresion(ruta_archivo='../data/Leukemia_GSE9476.csv'):
    try:
        df = pd.read_csv(ruta_archivo)
        return df
    except Exception as e:
        print(f"error al cargar los datos: {e}")
        return None

def analisis_sobreexpresion(df, id_gen, nombre_gen, archivo_salida):
    """
    genera un diagrama de caja para evaluar cómo varía la expresión
    de un gen específico a través de los diferentes subtipos de leucemia.
    """
    plt.figure(figsize=(10, 6))
    
    # el eje x será el tipo de leucemia, el eje y será el nivel de expresión del gen
    sns.boxplot(x='type', y=id_gen, data=df, boxprops={'alpha': 0.8})
    sns.stripplot(x='type', y=id_gen, data=df, color=".25", alpha=0.6, jitter=True)
    
    plt.title(f'perfil de expresión diferencial: {nombre_gen} ({id_gen})', fontsize=14)
    plt.xlabel('subtipo clínico de leucemia', fontsize=12)
    plt.ylabel('nivel de expresión (microarray)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(archivo_salida, dpi=300)
    plt.close()
    print(f"análisis de sobreexpresión guardado: {archivo_salida}")

def analisis_correlacion_cruzada(df, id_gen1, nombre_gen1, id_gen2, nombre_gen2, archivo_salida):
    """
    calcula el coeficiente de correlación de pearson y genera un gráfico de dispersión
    con línea de regresión para evaluar la co-expresión de dos genes.
    """
    datos_limpios = df[[id_gen1, id_gen2, 'type']].dropna()
    
    coeficiente, p_valor = pearsonr(datos_limpios[id_gen1], datos_limpios[id_gen2])
    
    plt.figure(figsize=(9, 7))
    
    # gráfico de dispersión con línea de regresión lineal
    sns.regplot(x=id_gen1, y=id_gen2, data=datos_limpios, scatter=False, color='gray')
    sns.scatterplot(x=id_gen1, y=id_gen2, hue='type', data=datos_limpios, s=60, alpha=0.8)
    
    plt.title(f'matriz de co-expresión: {nombre_gen1} vs {nombre_gen2}', fontsize=14)
    plt.xlabel(f'expresión de {nombre_gen1} ({id_gen1})', fontsize=12)
    plt.ylabel(f'expresión de {nombre_gen2} ({id_gen2})', fontsize=12)
    
    texto_estadistico = f"correlación de pearson (r): {coeficiente:.3f}\nvalor p: {p_valor:.2e}"
    plt.annotate(texto_estadistico, xy=(0.05, 0.95), xycoords='axes fraction', 
                 fontsize=11, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
                 verticalalignment='top')
    
    plt.legend(title='tipo de leucemia', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(archivo_salida, dpi=300)
    plt.close()
    
    print(f"análisis de correlación cruzada guardado: {archivo_salida}")
    print(f"estadística de red reguladora - pearson r: {coeficiente:.3f} (p-valor: {p_valor:.2e})")

if __name__ == "__main__":
    id_loc401317 = '205931_s_at'
    id_kcnj15 = '210119_at' # usamos el que obtuvo mayor f-score
    
    df_pacientes = cargar_datos_expresion()
    
    if df_pacientes is not None:
        print("\niniciando validación genómica...")
        
        # 1. análisis del lncrna loc401317
        analisis_sobreexpresion(
            df_pacientes, 
            id_loc401317, 
            "loc401317 (materia oscura / lncrna)", 
            "../plots/sobreexpresion_loc401317.png"
        )
        
        # 2. análisis del canal iónico kcnj15
        analisis_sobreexpresion(
            df_pacientes, 
            id_kcnj15, 
            "kcnj15 (canal de potasio / resistencia)", 
            "../plots/sobreexpresion_kcnj15.png"
        )
        
        # 3. análisis de correlación 
        analisis_correlacion_cruzada(
            df_pacientes,
            id_loc401317, "loc401317",
            id_kcnj15, "kcnj15",
            "../plots/correlacion_kcnj15_loc401317.png"
        )