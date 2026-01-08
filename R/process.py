import numpy as np
import pandas as pd
import formulaic
import itertools

def _one_indexed(df, colname):
    return df[colname].factorize()[0] + 1

def filter():
    pass

def py_process(df, control_treatment, G):
    df = df.copy()

    ## group RNA dispersion by mean count per batch
    df = _group_dispersion_indices(df, G)

    # create normalization factors
    df = _create_normalization_factors(df)

    # create effect model matrix
    df.loc[~df['is_ctrl'], 'architecture_model_matrix'] = df.loc[~df['is_ctrl'], 'architecture'] # exclude scramble/negative controls
    df.loc[df['is_ctrl'], 'architecture_model_matrix'] = np.nan
    df.loc[~(df['treatment'] == control_treatment), 'treatment_model_matrix'] = df.loc[~(df['treatment'] == control_treatment), 'treatment'] # exclude control treatment
    df.loc[(df['treatment'] == control_treatment), 'treatment_model_matrix'] = np.nan

    X, Xw, Xv, Xu, effects = _create_model_matrix_csr(df, 'architecture_model_matrix', 'treatment_model_matrix')

    # create correction factor model matrix
    if 'correction_factor' not in df.columns:
        df['correction_factor'] = 1
    Y, Yw, Yv, Yu, covariates = _create_model_matrix_csr(df, 'correction_factor', 'treatment_model_matrix')

    df['tre_id'] = _one_indexed(df, 'architecture')
    tres = df[['tre_id', 'architecture']].drop_duplicates().sort_values(by='tre_id')['tre_id'].tolist()

    data = {
        'P': len(correction_effects),
        'E': len(df['architecture'].unique()),
        'B': len(df['batch'].unique()),
        'N': df.shape[0],
        'N_EFFECTS': len(effects),
        'to_tre': df['tre_id'].tolist(),
        'R': df['RNA_count'].tolist(),
        'd': df['d_emb'].tolist(),
        'S_r': df['S_r'].tolist(),
        'to_dispersion_index': df['batch_dispersion_index'].tolist(),
        'N_DISPERSION_GROUPS': int(df['batch_dispersion_index'].max()),
        'n_w': len(Xw),
        'n_v': len(Xv),
        'n_u': len(Xu),
        'Xw': Xw.tolist(),
        'Xv': Xv.tolist(),
        'Xu': Xu.tolist(),
        'Xm': X.shape[0],
        'Xn': X.shape[1],
        'y_n_w': len(Yw),
        'y_n_v': len(Yv),
        'y_n_u': len(Yu),
        'Yw': Yw.tolist(),
        'Yv': Yv.tolist(),
        'Yu': Yu.tolist(),
        'Ym': Y.shape[0],
        'Yn': Y.shape[1],
    }

    keju = {'df': df, 'data': data, 'tres': tres, 'effects': effects, 'covariates': correction_effects}
    return keju
    
def py_use_motif_shrinkage(keju, motif='motif'):
    df = keju['df']
    data = keju['data']

    df['alpha_motif_id'] = _one_indexed(df, motif)
    tre_to_motif = df[['architecture', motif, 'tre_id', 'alpha_motif_id']].drop_duplicates()
    tre_to_motif = tre_to_motif.sort_values(by='tre_id')
    data['tre_to_motif'] = tre_to_motif['alpha_motif_id'].tolist()
    alpha_motifs = tre_to_motif[motif].tolist()

    # ignore negative control tres
    # df.loc[~df['is_ctrl'], 'beta_motif_id'] = _one_indexed(df.loc[~df['is_ctrl']], motif)
    beta_to_motif = df.loc[~df['is_ctrl'], ['architecture', motif]].drop_duplicates()
    motif_lookup = {
        x['architecture']: x[motif] for _, x in beta_to_motif.iterrows()
    }

    effects = data['effects']

    # create an arbitrary ordering/indexing 
    beta_motifs = beta_to_motif[motif].unique().tolist()
    data['effect_to_motif_x_treatment'] = [beta_motifs.index(motif_lookup[x]) + 1 for x in effects] # one indexed

    keju['alpha_motifs'] = alpha_motifs
    keju['beta_motifs'] = beta_motifs
    return keju

def py_use_covariate_slope_intercept(keju, motif='motif', covariate='covariate'):
    df = keju['df']
    data = keju['data']
    
    covariates = keju['covariates']
    df['covariate_id'] = df[covariate].map(lambda x: covariates.index(x) + 1)

    tre_to_covariate = df[['architecture', 'tre_id', covariate, 'covariate_id']].drop_duplicates()
    tre_to_covariate = tre_to_covariate.sort_values(by='tre_id')
    data['tre_to_correction_factor'] = tre_to_covariate['covariate_id'].tolist()
    keju['covariates'] = covariates
    return keju

def _create_normalization_factors(df):
    S_d = df.groupby('dna_batch')['DNA_count'].quantile(.75)
    df['d_emb'] = df.apply(lambda x: x['DNA_count'] / S_d[x['dna_batch']], axis=1)
    df['S_d'] = df['dna_batch'].map(lambda x: S_d[x])

    S_r = df.groupby('batch')['RNA_count'].quantile(.75)
    df['S_r'] = df['batch'].apply(lambda x: S_r[x])
    return df

def _group_dispersion_indices(df, G, factors = ['batch']):
    factor_instances = {f: list(df[f].unique()) for f in factors}
    factor_products = list(itertools.product(*factor_instances.values()))
    
    rna_sorted = df.groupby(factors + ['architecture'])['RNA_count'].mean().sort_values().reset_index()
    rna_sorted = rna_sorted.sort_values(by=factors + ['RNA_count'])
    rna_sorted = rna_sorted.set_index(factors)

    i = 0

    for f in df[factors].drop_duplicates().iterrows():
        mini_factors = tuple(f[1].values)

        rna_sorted.loc[mini_factors, 'dispersion_index'] = range(i, i + rna_sorted.loc[mini_factors].shape[0])                
        rna_sorted.loc[mini_factors, 'dispersion_index'] = (rna_sorted.loc[mini_factors, 'dispersion_index'] / G).astype(int) + 1
        i = int(G * rna_sorted['dispersion_index'].max())

    rna_sorted['dispersion_index'] = rna_sorted['dispersion_index'].astype(int)

    dispersion_mapping = {
        (batch, row['architecture']): row['dispersion_index'] for batch, row in rna_sorted.iterrows()
    }

    df['batch_dispersion_index'] = df.apply(lambda x: dispersion_mapping[(x.batch, x.architecture)], axis=1).astype(int)

    return df

def _create_model_matrix_csr(df, col1, col2):
    X = formulaic.Formula(f'~0 + {col1} : {col2}').get_model_matrix(
                                                                    df, 
                                                                    ensure_full_rank=False,
                                                                    na_action='ignore',
                                                                    output='sparse'
                    )
    X_csr = X.tocsr()

    Xw = X_csr.data
    Xv = X_csr.indices + 1
    Xu = X_csr.indptr + 1

    effects = X.model_spec.column_names
    effects = [x.split(f']:{col2}[')[0].split(f'{col1}[')[1] for x in effects]
    # effects = [x.split(f']:{col2}[T.')[0].split(f'{col1}[T.')[1] for x in effects] # reticulate throws error because of T. contrast levels

    return X, Xw, Xv, Xu, effects