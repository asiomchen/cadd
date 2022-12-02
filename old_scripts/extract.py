import os
import re

import pandas as pd


def extract(filename):
    '''
    Extracts information from the docking results
    :param filename: smina or gnina output file .sdf
    :return: pd.Series with the extracted information
    '''
    # extract line started with ">  <minimizedAffinity>" and the next line
    compound_db_name_pattern = re.compile(r'CHEMBL\d+')
    compound_db_name = compound_db_name_pattern.search(filename)
    if compound_db_name is None:
        compound_db_name = filename.split('.')[0]
    else:
        compound_db_name = compound_db_name.group()

    result = []
    with open(filename) as f:
        for line in f:
            if re.match(r">  <minimizedAffinity>", line):
                result.append(next(f))
            if re.match(r"> <minimizedAffinity>", line):
                result.append(next(f))
            if re.match(r">  <CNNscore>", line):
                result.append(next(f))
            if re.match(r">  <CNNaffinity>", line):
                result.append(next(f))
    # delete the last newline character
    result = [x.strip() for x in result]
    # convert to float
    result = [float(x) for x in result]
    result.insert(0, compound_db_name)
    return pd.Series(result)


if __name__ == '__main__':
    os.chdir('/home/anton/in_dev/Docking_tools/docked_smina_auto')
    files = os.listdir()
    df = pd.DataFrame()
    for file in files:
        df = df.append(extract(file), ignore_index=True)
    # df.columns = ['ID', 'Affinity_1', 'CNNscore_1', 'CNNaffinity_1', 'Affinity_2', 'CNNscore_2', 'CNNaffinity_2', 'Affinity_3', 'CNNscore_3', 'CNNaffinity_3', 'Affinity_4', 'CNNscore_4', 'CNNaffinity_4', 'Affinity_5', 'CNNscore_5', 'CNNaffinity_5']
    # df.sort_values(by=['Affinity_1'], inplace=True, ascending=True)
    # df = df[['ID', 'Affinity_1', 'CNNscore_1', 'CNNaffinity_1']]
    # df['Realative_affinity'] = df['Affinity_1'] * df['CNNscore_1']
    # df.sort_values(by=['Realative_affinity'], inplace=True, ascending=True)
    # print(df)
    #
    # # smina part
    smina_columns = ['ID', 'Affinity_1', 'Affinity_2', 'Affinity_3', 'Affinity_4', 'Affinity_5']
    df.columns = smina_columns
    df.sort_values(by=['Affinity_1'], inplace=True, ascending=True)
    print(df)
    df.to_csv('/home/anton/in_dev/Docking_tools/smina_128_results.csv', index=False)
