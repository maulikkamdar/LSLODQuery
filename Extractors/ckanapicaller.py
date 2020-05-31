#!/usr/bin/python
# -*- coding: utf-8 -*-
# CKAN Calls -- not required

import networkx as nx
from utils import HTTPUtils, MatrixIO, FileUtils
import urllib

data_hub_url = 'https://old.datahub.io/api/3/action/'


def is_ascii(string):
    return all(ord(char) < 128 for char in string)


def query_ckan(url):
    hu = HTTPUtils()
    response = hu.get_json(url, None)  # No API Key required
    return response


def get_all_tag_list():
    url = data_hub_url + 'tag_list'
    tag_list = query_ckan(url)
    ls_tag_list = []
    all_tag_s = ['bio', 'life', 'medic', 'clinic']
    for k in tag_list['result']:

        # process on terms with life or bio in them ...

        for m in all_tag_s:
            if m in k:
                if not is_ascii(k):
                    continue
                ls_tag_list.append(k)
    return ls_tag_list


def ret_ckan_datasets(tag):
    url = data_hub_url + 'package_search?fq=tags:' + urllib.quote(tag) \
        + '&rows=500'  # lifesciences
    print url
    response = query_ckan(url)

    # dataset_list = response
    # return dataset_list

    dataset_list = {}
    for k in response['result']['results']:
        if k['name'] not in dataset_list:
            dataset_list[k['name']] = k
    return dataset_list


def ret_ckan_dataset_info(dataset):

    # query for a specific dataset url\

    return None


ls_tag_list = get_all_tag_list()

# lifescience_dset_list = ret_ckan_datasets("lifesciences")

int_dset = {}
for tag in ls_tag_list:
    _dset_list = ret_ckan_datasets(tag)
    print 'returning datasets for ' + tag + ' : ' + str(len(_dset_list))
    for m in _dset_list:
        if m not in int_dset:
            int_dset[m] = _dset_list[m]

print int_dset['bio2rdf-drugbank']

ref_set = {}
for k in int_dset:
    res_set = []
    for m in int_dset[k]['resources']:
        res_set.append({'format': m['format'], 'url': m['url'],
                       'name': m['name']})
    tag_set = []
    for m in int_dset[k]['tags']:
        tag_set.append(m['name'])
    ref_set[k] = {
        'tags': tag_set,
        'resources': res_set,
        'title': int_dset[k]['title'],
        'url': int_dset[k]['url'],
        'extras': int_dset[k]['extras'],
        'notes': int_dset[k]['notes'],
        }

lslodG = nx.DiGraph()
for k in ref_set:
    lslodG.add_node(k)
    for m in ref_set[k]['extras']:
        key = m['key'].split(':')
        if key[0] == 'links':
            lslodG.add_edge(k, key[1], weight=m['value'])

nx.write_graphml(lslodG, 'lslodG_ckan.graphml')
mfio = MatrixIO()
mfio.save_matrix(ref_set, 'lslod_dset_ckan.dat')
for k in ref_set:
    spr = ''
    for m in ref_set[k]['resources']:
        if m['format'] == 'api/sparql' or 'sparql' in m['url']:
            spr = m['url']
    print (k, ref_set[k]['title'], spr)
