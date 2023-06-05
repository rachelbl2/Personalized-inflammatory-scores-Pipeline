import pandas as pd
from statsmodels import api as sm
from sklearn.preprocessing import StandardScaler
import numpy as np

def expression_preproccesing(expression_df: pd.DataFrame,control_individuals:list):
    '''
    :param full_data: full expression data set. index - genes, columns - individuals
    :param control_individuals: list of columns used for reference data
    :return: the normalized expression matrix
    '''
    expression_df = np.log2(expression_df+1)
    return normalize_data(expression_df[expression_df[control_individuals]],expression_df) #normalize based on contol individuals

def normalize_data(ref_data:pd.DataFrame, data:pd.DataFrame)->pd.DataFrame:
    '''
    noramalize df based on reference data
    :param data: the full data frame. index are genes and columns individuals
    :param ref_data: reference data frame to normalize on - a subset of the full data.
    index are genes and columns individuals
    :return: the normalized full data frame
    '''
    scaler1 = StandardScaler()
    scaler2 = StandardScaler()

    data_scaled = pd.DataFrame(scaler1.fit_transform(data), index=data.index, columns=data.columns)
    ref_data_scaled = pd.DataFrame(scaler1.fit_transform(ref_data),
                                   index=ref_data.index, columns=ref_data.columns)
    scaler2.fit(ref_data_scaled.T)
    data_scaled = pd.DataFrame(scaler2.transform(data_scaled.T).T,
                               index=data_scaled.index, columns=data_scaled.columns)
    return data_scaled


def calculate_state(df:pd.DataFrame,x:str,y:str,individual:str,const: bool = True):
    '''
    calculates state of certain individual
    :param df: data frame with map position columns and the individual's expression column
    :param x: column name for map x positions
    :param y: column name for map y positions
    :param individual: column name of the individual's gene expression
    :param const: if True adds const to the regression
    '''
    if const:
        X = sm.add_constant(df[[x, y]])
    else:
        X = df[[x, y]]
    res = sm.OLS(df[individual], X).fit()
    return res, res.params[x], res.params[y]

def center_map(map:pd.DataFrame,x:str,y:str)->pd.DataFrame:
    '''
    center map positions around 0 (median)
    :param map: df with the map x,y positions
    :param x: column name for map x positions
    :param y: column name for map y positions
    :return: map centered around 0
    '''
    median_x = map[x].median()
    median_y = map[y].median()
    centered_map = map.copy()
    centered_map[x] = centered_map[x] - median_x
    centered_map[y] = centered_map[y] - median_y
    return centered_map


def states_df(map:pd.DataFrame,x:str,y:str,df:pd.DataFrame,const:bool = True,center :bool= True)->pd.DataFrame:
    '''
    calculates state for every individual in df based on map with (x,y) axis.
    :param map: gene space positions (T/R  or MetS/SI)
    :param df: normalized data frame in which the columns are individuals and the index are genes
    :param const: if True adds const to the regression
    :param center: if True map coordinates are centered around 0
    :return: Dataframe with the S1 (horizontal axis) and S2 (vertical axis) states and the corresponding p values
    '''
    if center:
        map = center_map(map,x,y)
    states_df = pd.DataFrame(columns=['individual',f'{x} state',f'{y} state',f'pval {x}',f'pval {x}'])
    individuals = df.columns
    df = df.join(map).dropna()
    for individual in individuals:
        res, res_x, res_y= calculate_state(df,x,y,individual,const)
        pvals = res.pvalues
        states_df = states_df.append({'individual':individual, f'{x} state': res_x, f'{y} state': res_y,
                                      f'pval {x}':pvals[0],f'pval {y}':pvals[1]},ignore_index=True)

    states_df = states_df.set_index('individual')
    return states_df


def norm_state_on_control(states_df: pd.DataFrame,state_name: str,contol_individuals: list) -> pd.DataFrame:
    '''
    standardization of the states based on the controls - subtracting the mean level (centering) and
    dividing by the standard deviation of the control subjects
    :param state_name: the name of the state column to normalize
    :param contol_individuals: list of control individuals
    :return: the state df with the state_name column standardize based on the controls subjects
    '''
    controls_states = states_df[states_df.index.isin(contol_individuals)]
    mean = controls_states[state_name].mean()
    std = controls_states[state_name].std()
    norm_states_df = states_df.copy()
    norm_states_df[state_name] = (norm_states_df[state_name] - mean)/std
    return norm_states_df

def norm_all_states(states_df: pd.DataFrame,control_individuals: list)-> pd.DataFrame:
    '''
    standardization of the states based on the controls for all the states columns (R, SI, T and MetS)
    '''
    states_df = norm_state_on_control(states_df,'SI',control_individuals)
    states_df = norm_state_on_control(states_df,'R',control_individuals)
    states_df = norm_state_on_control(states_df,'T',control_individuals)
    states_df = norm_state_on_control(states_df,'Mets',control_individuals)
    return states_df

def R_SI_balance_score(full_states_table: pd.DataFrame) -> pd.DataFrame:
    '''
    calculates the R/SI balance score by subtracting SI from R (can perform only after the
    standarization of R and SI levels).
    :param full_states_table: normalized states table
    :return: states table with R/SI balance score column
    '''
    full_states_table[f'R/SI balance score'] = full_states_table.apply(lambda row: row['R'] - row['SI'],axis=1) #use only after the normalization of the states (norm_all_states)
    return full_states_table

def full_R_SI_score_pipeline(expression_df,control_individuals):
    '''
    runs the entire pipline from gene expression matrix to full states table with R/SI balance score column
    :param expression_df: gene expression table in which the columns are individuals and the index are genes
    :param control_individuals: list of controls individuals (subset of columns)
    '''
    T_R_map = pd.read_csv('resistance_tolerance_map.csv',index_col=0)
    mets_SI_map = pd.read_csv('MetS_systemic_inflammation_map.csv',index_col=0)
    expression_df = expression_preproccesing(expression_df,control_individuals)
    r_t_states = states_df(T_R_map,'T','R',expression_df)
    mets_si_states = states_df(mets_SI_map,'MetS','SI',expression_df)
    full_states_table = r_t_states.join(mets_si_states)
    full_states_table = norm_all_states(full_states_table,control_individuals)
    full_states_table = R_SI_balance_score(full_states_table)
    return full_states_table
