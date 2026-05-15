import os
import pandas as pd
import numpy as np
from scipy.stats import zscore

os.environ.pop('mplbackend', None)
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

def generar_mapa_calor(ruta_expresion, ruta_biomarcadores, ruta_diccionario, archivo_salida):
    """
    genera un heatmap estandarizado cruzando los top biomarcadores con los subtipos de leucemia.
    """
    try:
        df_pacientes = pd.read_csv(ruta_expresion)
        df_biomarcadores = pd.read_csv(ruta_biomarcadores)
        df_dict = pd.read_csv(ruta_diccionario, sep='\t', comment='#', low_memory=False)
        
        top_ids = df_biomarcadores['gen_id'].head(20).tolist()
        
        columnas_filtro = ['type'] + top_ids
        df_filtrado = df_pacientes[columnas_filtro]
        
        # calcular la media de expresión por subtipo de leucemia
        df_medias = df_filtrado.groupby('type').mean()
        
        # mapear nombres reales a los genes
        nombres_columnas = []
        for gen_id in df_medias.columns:
            match = df_dict[df_dict['ID'] == gen_id]['Gene Symbol'].values
            nombre = str(match[0]) if len(match) > 0 and pd.notna(match[0]) else f"materia_oscura ({gen_id})"
            nombres_columnas.append(nombre)
            
        df_medias.columns = nombres_columnas
        
        # estandarización (z-score)
        # esto es vital: compara genes de alta expresión con genes de baja expresión en la misma escala
        df_zscore = df_medias.apply(zscore, axis=0)
        
        plt.figure(figsize=(14, 8))
        sns.heatmap(df_zscore, cmap="coolwarm", center=0, annot=True, fmt=".1f", 
                    linewidths=.5, cbar_kws={'label': 'expresión relativa (z-score)'})
        
        plt.title('perfil de expresión molecular por subtipo de leucemia', fontsize=14, pad=20)
        plt.xlabel('biomarcadores (top 20)', fontsize=12)
        plt.ylabel('subtipo clínico', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        plt.savefig(archivo_salida, dpi=300)
        plt.close()
        print(f"mapa de calor generado con éxito en: {archivo_salida}")

    except Exception as e:
        print(f"error durante la generación del mapa de calor: {e}")

if __name__ == "__main__":
    ruta_datos = '../data/Leukemia_GSE9476.csv'
    ruta_biomarcadores = '../results/biomarcadores.csv'
    ruta_anotacion = '../data/GPL570-55999.txt'
    ruta_salida = '../plots/heatmap_subtipos_top20.png'
    
    generar_mapa_calor(ruta_datos, ruta_biomarcadores, ruta_anotacion, ruta_salida)