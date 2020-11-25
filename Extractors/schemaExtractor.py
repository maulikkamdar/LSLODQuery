#!/usr/bin/python
# -*- coding: utf-8 -*-
import json
import re
import sys
import os
import subprocess
import time
import urllib
import rdflib
import hashlib
from operator import itemgetter
from SPARQLWrapper import SPARQLWrapper, JSON, POSTDIRECTLY
import networkx as nx
from utils import MatrixIO
import xml.etree.ElementTree as ET
import argparse

first_cap_re = re.compile('(.)([A-Z][a-z]+)')
all_cap_re = re.compile('([a-z0-9])([A-Z])')


def get_endpoint(endpoint):
    g = rdflib.Graph('TPFStore')
    g.open(endpoint)
    return g


def get_simple_results(endpoint, query):
    g = get_endpoint(endpoint)
    results = g.query(query)
    return results


def get_simple_results_ldf(endpoint, query):
    proc = subprocess.Popen(['/Users/mkamdar/Client.js/bin/ldf-client',
                             endpoint, query], stdout=subprocess.PIPE)
    output = proc.stdout.read()
    return output


def proc_floats(val):
    val_parts = val.split('^^')
    if len(val_parts) > 1:
        return float((val_parts[0])[1:-1])
    else:
        return 0


def get_simple_results_normal(endpoint, query, limres=False):
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    sparql.setTimeout(3600000)

    # sparql.setCredentials("admin", "admin")
    # if limres: sparql.addExtraURITag("limit", "100")

    try:
        results = sparql.query().convert()
    except Exception:

        # print results
        # results = process_xml_results(results)

        print('Error encountered for ' + query)
        return None
    return results


def parse_camel_case(name):
    s1 = first_cap_re.sub(r'\1_\2', name)
    return all_cap_re.sub(r'\1_\2', s1).lower()


def parse_uri_token(token):
    prop = parse_camel_case(token)
    pparts = re.split('[-_]', prop)
    name = ' '.join([x.title() for x in pparts])
    return name


def parseURI(uri):
    nparts = re.split('[#:/]', uri)
    name = parse_uri_token(nparts[len(nparts) - 1])
    return name


def classlist_refined(qe, inc_counts=True):
    if inc_counts:
        q1 = \
            'SELECT ?Concept (COUNT (?x) AS ?cCount) WHERE {?x a ?Concept} GROUP BY ?Concept ORDER BY DESC(?cCount)'
    else:
        # if you run into a timeout error
        q1 = 'SELECT DISTINCT * WHERE {[] a ?Concept}'
    res1 = get_simple_results_normal(qe, q1)
    class_list = {}
    if not res1:
        if inc_counts:
            class_list = classlist_refined(qe, inc_counts=False)
            return class_list
        else:
            return None  # to implement revert to other query
    for result in res1['results']['bindings']:
        con_uri = str(result['Concept']['value'].encode('utf-8'))
        if con_uri[0:4] != 'http':  # because non-standard uris
            continue
        if 'openlinksw' in con_uri:  # because virtuosos
            continue
        print(con_uri)
        class_count = (int(result['cCount']['value'
                                            ]) if inc_counts else 0)
        class_list[con_uri] = class_count
    print(len(class_list))
    return class_list


def process_xml_results(results):
    a = results.toxml().replace('\n', '')
    axml = ET.fromstring(a)
    json_results = {}
    json_results['results'] = {'bindings': []}
    for res in \
            axml.iter('{http://www.w3.org/2005/sparql-results#}result'):
        binding = {}
        for child in res:
            binding[child.attrib['name']] = {'value': child[0].text}
        json_results['results']['bindings'].append(binding)
    return json_results


def property_refined(qe, class_list, version=1):

    # assume to get all properties first, save unfound classes and rerun on a reduced randomized set and estimate instances

    G = nx.DiGraph()
    under_query = []
    if version == 1:
        q3_template = \
            "SELECT DISTINCT ?p ?c (COUNT(?x) AS ?count) ?valType WHERE {?x a <CONCEPT_URI>; ?p ?o . OPTIONAL {?o a ?c} . FILTER(!(?p = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type')) . BIND(DATATYPE(?o) AS ?valType)} GROUP BY ?p ?c ?valType ORDER BY DESC(?count)"
    elif version == 3:
        q3_template = \
            "SELECT DISTINCT ?p ?c (COUNT(?x) AS ?count) ?valType WHERE {?x a <CONCEPT_URI>; ?p ?o . OPTIONAL {?o a ?c} . FILTER(!(?p = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type')) . BIND(DATATYPE(?o) AS ?valType)} ORDER BY DESC(?count)"
    else:
        q3_template = \
            "SELECT DISTINCT ?p ?c (COUNT(?x) AS ?count) WHERE {?x a <CONCEPT_URI>; ?p ?o . OPTIONAL {?o a ?c} . FILTER(!(?p = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type')) } GROUP BY ?p ?c ORDER BY DESC(?count)"  # BIND is not supported!
    for k in class_list:
        avoid_namespace = set(['purl.obolibrary.org',
                               'ncicb.nci.nih.gov',
                               'purl.bioontology.org'])
        kparts = k.split('/')
        if len(kparts) > 2:
            print((k, kparts[2], class_list[k]))
            if kparts[2] in avoid_namespace and class_list[k] < 2:
                continue
        time.sleep(5)
        G.add_node(k, type='class')
        count = 0
        q3 = q3_template.replace('CONCEPT_URI', k)
        print(q3)
        res3 = get_simple_results_normal(qe, q3, limres=True)
        if not res3:
            under_query.append(k)
            continue
        if isinstance(res3, basestring):
            under_query.append(k)
            continue
        for inst_res in res3['results']['bindings']:
            pred = inst_res['p']['value']
            targ = (inst_res['c']['value'] if 'c' in inst_res else None)
            count = int(inst_res['count']['value'])
            valType = (inst_res['valType']['value'] if 'valType'
                       in inst_res else '')
            if targ:
                G.add_node(targ, type='class')
                if not G.has_edge(k, targ):
                    G.add_edge(k, targ, type=[pred], count=[count])
                else:
                    G[k][targ]['type'].append(pred)
                    G[k][targ]['count'].append(count)
            else:
                G.add_node(pred, type='literal')
                G.add_edge(k, pred, type=valType, count=count)
        nx.write_gpickle(G, map_name + '.gpickle')
    return (G, under_query)


def init_main(args):
    qe = args.qe
    map_name = args.name
    print(("Executing for", qe, map_name))
    mfio = MatrixIO()
    _classList = classlist_refined(qe)
    mfio.save_matrix(_classList, '_classList_' + map_name + '.dat')

    (G, under_query) = property_refined(qe, _classList, version=1)
    nx.write_gpickle(G, map_name + '.gpickle')
    mfio.save_matrix(under_query, 'under_query_' + map_name + '.dat')

    for k in _classList:
        G.node[k]['count'] = _classList[k]

    nx.write_gpickle(G, map_name + '.gpickle')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        # don't want to reformat the description message, use file's doc_string
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    parser.add_argument("--qe", help="Query Endpoint")
    parser.add_argument("--name", help="Map Name")
    args = parser.parse_args()
    init_main(args)
