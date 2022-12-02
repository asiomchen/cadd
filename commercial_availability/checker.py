import sys
import warnings

import requests

sys.path.insert(1, './')


class Molport:
    def __init__(self, smiles=None, molport_id=None, is_commercial=None):
        self.smiles = smiles
        self.is_commercial = is_commercial
        self.molport_id = molport_id
        self.username = "john.spade"
        self.password = "fasdga34a3"
        self.payload = {
            "User Name": self.username,
            "Authentication Code": self.password,
            "Structure": self.smiles,
            "Search Type": "EXACT",
            "Maximum Search Time": 60000,
            "Maximum Result Count": 1000,
            "Chemical Similarity Index": 0.9
        }

    def find_compound(self, smiles=None):
        """
        Finds the Molport ID of a compound. If compound have molport ID exists,
         assupms that it is commercial.
        :param smiles: canonical smiles string
        :return:
        """
        payload = {
           "User Name": self.username,
           "Authentication Code": self.password,
           "Structure": smiles,
           "Search Type": 4,
           "Maximum Search Time": 60000,
           "Maximum Result Count": 10000,
           "Chemical Similarity Index": 1
        }
        similarity_request = requests.post('https://api.molport.com/api/chemical-search/search', json=payload)
        response = similarity_request.json()
        try:
            self.molport_id = response['Data']['Molecules'][0]['MolPort Id'][8:]
            self.is_commercial = True
        except:
            self.molport_id = None
            self.is_commercial = False
        return Molport(smiles=smiles, molport_id=self.molport_id, is_commercial=self.is_commercial)


    def get_compound_suppliers(self, smiles=None, as_df=False):
        self.smiles = smiles
        if self.molport_id is None:
            self.molport_id = self.find_compound(self.smiles)
        if not self.is_commercial:
            warnings.warn("This compound is non-commercial and cannot be retrieved from Molport")
            return None
        molport_id_request = 'https://api.molport.com/api/molecule/load?' \
                             'molecule={}' \
                             '&username=john.spade' \
                             '&authenticationcode=fasdga34a3'

        r2 = requests.get(molport_id_request.format(self.molport_id))
        response = r2.json()
        results = response['Data']['Molecule']['Catalogues']['Screening Block Suppliers']
        if as_df:
            df = pd.DataFrame()
            for supplier in results:
                df = df.append(supplier, ignore_index=True)
            shipping_options = pd.DataFrame()
            for s_cost, supplier in zip(df['Shipment Costs'], df['Supplier Name']):
                shipping_option = pd.DataFrame(s_cost, index=[supplier for i in range(len(s_cost))])
                shipping_options = shipping_options.append(shipping_option)
            catalogs = pd.DataFrame()
            for s_cost, supplier in zip(df['Catalogues'], df['Supplier Name']):
                catalog = pd.DataFrame(s_cost, index=[supplier for i in range(len(s_cost))])
                catalogs = catalogs.append(catalog)

            merged = pd.merge(shipping_options, catalogs, left_index=True, right_index=True, )
            return merged
        else:
            return results



# searcher = Molport()
# molport_id = searcher.find_compound(smiles="NCC1OC(OC2C(N)CC(N)C(O)C2OC2OC(CO)C(OC3OC(CN)C(O)C(O)C3O)C2O)C(N)C(O)C1O")
# print(searcher.is_commercial)
# results = searcher.get_compound_suppliers(
#     smiles="NCC1OC(OC2C(N)CC(N)C(O)C2OC2OC(CO)C(OC3OC(CN)C(O)C(O)C3O)C2O)C(N)C(O)C1O", as_df=True)
# print(results)
# results.to_csv("molport_results.csv")

if __name__ == '__main__':
    search_unavailable = Molport().find_compound(smiles='CCCC').is_commercial  # False
    search_available = Molport().find_compound(smiles='CC(=O)Nc1ccc(cc1)O').is_commercial  # True
    import pandas as pd
    df = pd.read_csv('/home/anton/in_dev/Docking_tools/master/cox_chembl_cleaned_v2.csv')
    df['is_commercial'] = df['Smiles'].apply(lambda x: Molport().find_compound(smiles=x).is_commercial)
    print(df['is_commercial'].value_counts())
    df.to_csv('cox_chembl_cleaned_v2.csv')


