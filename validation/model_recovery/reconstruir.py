"""Reconstruye artefactos sin entrenar ni modificar el proyecto original.
Uso: python reconstruir.py --proyecto RUTA --salida CARPETA
Dependencias: numpy, pandas, scipy, scikit-learn, h5py, joblib.
"""
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import h5py
import joblib
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder, StandardScaler


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--proyecto', required=True, type=Path)
    ap.add_argument('--salida', required=True, type=Path)
    args = ap.parse_args()
    root, out = args.proyecto.resolve(), args.salida.resolve()
    if out == root or root in out.parents:
        raise ValueError('La salida debe estar fuera del proyecto original')
    out.mkdir(parents=True, exist_ok=True)
    matrix = root / 'data/Leukemia_GSE9476.csv'
    model = root / 'models/modelo_leucemia.keras'
    biomarkers = root / 'results/biomarcadores.csv'
    df = pd.read_csv(matrix)
    X = df.drop(columns=['samples', 'type'])
    encoder = LabelEncoder().fit(df.type)
    y = encoder.transform(df.type)
    b = pd.read_csv(biomarkers)
    comparisons, splits, matches = [], [], []
    for fold, (tr, te) in enumerate(StratifiedKFold(5, shuffle=True, random_state=42).split(X, y), 1):
        scaler = StandardScaler().fit(X.iloc[tr])
        selector = SelectKBest(f_classif, k=500).fit(scaler.transform(X.iloc[tr]), y[tr])
        top = np.argsort(-selector.scores_, kind='stable')[:len(b)]
        same = X.columns[top].tolist() == b.gen_id.tolist()
        error = float(np.max(np.abs(selector.scores_[X.columns.get_indexer(b.gen_id)] - b.f_score.to_numpy())))
        comparisons.append(dict(fold=fold, top20_exact=same, max_f_score_error=error))
        splits.append(dict(fold=fold, train_row_indices=tr.tolist(), test_row_indices=te.tolist(),
                           train_samples=df.samples.iloc[tr].tolist(), test_samples=df.samples.iloc[te].tolist()))
        if same and error < 1e-8:
            matches.append((fold, tr, te, scaler, selector))
    if len(matches) != 1:
        raise RuntimeError(f'No hay coincidencia unica: {comparisons}')
    fold, tr, te, scaler, selector = matches[0]
    weights = {}
    with zipfile.ZipFile(model) as z:
        metadata = json.loads(z.read('metadata.json'))
        with h5py.File(io.BytesIO(z.read('model.weights.h5')), 'r') as h:
            for i, layer in enumerate(['dense', 'dense_1', 'dense_2']):
                weights[f'W{i}'] = h[f'layers/{layer}/vars/0'][()]
                weights[f'b{i}'] = h[f'layers/{layer}/vars/1'][()]
            optimizer_iterations = int(h['optimizer/vars/0'][()])
    selected = X.columns[selector.get_support()].to_numpy(dtype=str)
    # Dense/ReLU/Dropout/Dense/ReLU/Dropout/Dense/Softmax; dropout desactivado en inferencia.
    a = selector.transform(scaler.transform(X)).astype(np.float32)
    for i in range(3):
        a = a @ weights[f'W{i}'] + weights[f'b{i}']
        if i < 2:
            a = np.maximum(a, 0)
    a -= a.max(axis=1, keepdims=True)
    probs = np.exp(a)
    probs /= probs.sum(axis=1, keepdims=True)
    pred = probs.argmax(axis=1)
    predictions = pd.DataFrame(dict(sample=df.samples, real=df.type,
                                   predicted=encoder.inverse_transform(pred),
                                   partition=['test' if i in te else 'train' for i in range(len(df))]))
    for i, label in enumerate(encoder.classes_):
        predictions['prob_' + label] = probs[:, i]
    predictions.to_csv(out / 'predicciones_candidato.csv', index=False)
    pd.DataFrame({'position_0_based': np.arange(500), 'probe': selected}).to_csv(out / 'sondas_500_ordenadas.csv', index=False)
    joblib.dump(dict(scaler=scaler, selector=selector, encoder=encoder, input_columns=X.columns.tolist()), out / 'preprocesamiento_reconstruido.joblib')
    np.savez_compressed(out / 'modelo_inferencia_numpy.npz', **weights, mean=scaler.mean_, scale=scaler.scale_,
                        selected_indices=selector.get_support(indices=True), input_columns=X.columns.to_numpy(dtype=str),
                        selected_probes=selected, classes=encoder.classes_.astype(str))
    report = dict(status='preprocesamiento candidato reconstruido; no reentrenado',
                  matched_fold=fold, comparisons=comparisons, splits=splits,
                  classes=encoder.classes_.tolist(), model_metadata=metadata, optimizer_iterations=optimizer_iterations,
                  test_n=len(te), test_correct=int((pred[te] == y[te]).sum()),
                  candidate_test_accuracy=float(accuracy_score(y[te], pred[te])),
                  candidate_test_confusion_matrix=confusion_matrix(y[te], pred[te], labels=np.arange(len(encoder.classes_))).tolist(),
                  caveats=['La coincidencia del ranking identifica el fold del selector; no demuestra por si sola la procedencia de los pesos.',
                           'El resultado de un fold no verifica el 98.3% de validacion cruzada.',
                           'Inferencia calculada con NumPy a partir de los pesos; no comparada con TensorFlow.',
                           'No se han generado modelos de los otros cuatro folds ni reentrenado la red.'],
                  source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [matrix, model, biomarkers]})
    (out / 'verificacion.json').write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['matched_fold','test_n','test_correct','candidate_test_accuracy','optimizer_iterations']}, indent=2))


if __name__ == '__main__':
    main()
