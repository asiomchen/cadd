import os, sys
import subprocess
import urllib
import requests


'''
GET: https://data.plgrid.pl/list/[site_key]/[folder_path]
'''
DEFAULT_SITE_KEY = 'prometheus'
path = 'net/scratch/people/plgasiomchen/smina_test/Molecular_docking/results.tar.gz'
def get_data(site_key, folder_path):
    cmd = 'curl -X GET https://data.plgrid.pl/download/{}/{} --data-urlencode proxy="`cat grid_proxy`" -o downloaded.tar.gz'.format(site_key, folder_path, )
    print(cmd)

if __name__ == '__main__':
    get_data(DEFAULT_SITE_KEY, path)
