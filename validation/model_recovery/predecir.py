"""Inferencia con los pesos originales y el preprocesamiento candidato.
Uso: python predecir.py --csv matriz.csv --salida predicciones.csv
Requiere numpy y pandas. No valida utilidad clinica ni rendimiento externo.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def predict(frame, artifact):
    columns = artifact['input_columns'].tolist()
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f'Faltan {len(missing)} sondas; primeras: {missing[:10]}')
    x = frame[columns].to_numpy(dtype=np.float64)
    if not np.isfinite(x).all():
        raise ValueError('La expresion contiene valores ausentes o no finitos')
    x = ((x - artifact['mean']) / artifact['scale'])[:, artifact['selected_indices']].astype(np.float32)
    for i in range(3):
        x = x @ artifact[f'W{i}'] + artifact[f'b{i}']
        if i < 2:
            x = np.maximum(x, 0)
    x -= x.max(axis=1, keepdims=True)
    p = np.exp(x)
    p /= p.sum(axis=1, keepdims=True)
    result = pd.DataFrame({'sample': frame['samples'] if 'samples' in frame else np.arange(len(frame)),
                           'predicted': artifact['classes'][p.argmax(axis=1)]})
    for i, label in enumerate(artifact['classes']):
        result['prob_' + label] = p[:, i]
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True, type=Path)
    ap.add_argument('--modelo', type=Path, default=Path(__file__).parent / 'modelo_inferencia_numpy.npz')
    ap.add_argument('--salida', required=True, type=Path)
    args = ap.parse_args()
    if args.salida.resolve() in (args.csv.resolve(), args.modelo.resolve()):
        raise ValueError('La salida no puede sobrescribir una entrada')
    with np.load(args.modelo, allow_pickle=False) as artifact:
        result = predict(pd.read_csv(args.csv), artifact)
    result.to_csv(args.salida, index=False)
