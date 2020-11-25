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
from operator import itemgetter
from SPARQLWrapper import SPARQLWrapper, JSON
import networkx as nx
from copy import deepcopy
from utils import MatrixIO, FileUtils
import xml.etree.ElementTree as ET


def get_simple_results_normal(endpoint, query):
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    sparql.addExtraURITag('limit', '2000')
    #sparql.setCredentials('admin', 'admin')
    try:
        results = sparql.query().convert()
    except Exception, inst:
        print 'Error encountered for ' + query
        return None
    return results


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


def ret_inst(endpoint, node):
    query = 'SELECT ?x WHERE {?x a <' + node \
        + '>} ORDER BY RAND() LIMIT 2000'

    results = get_simple_results_normal(endpoint, query)
    if not str(type(results)) == "<type 'dict'>":
        return []
    if not results:
        return []
    instances = []
    for result in results['results']['bindings']:
        inst = result['x']['value']
        instances.append(inst)
    return instances


def ret_rels(endpoint, cnode, pnode):
    query = 'SELECT ?x WHERE {?c a <' + cnode + '>; <' + pnode \
        + '> ?x} ORDER BY RAND() LIMIT 2000'

    results = get_simple_results_normal(endpoint, query)
    if not str(type(results)) == "<type 'dict'>":
        return []
    if not results:
        return []
    instances = []
    for result in results['results']['bindings']:
        inst = result['x']['value']
        instances.append(inst)
    return instances


def query_instances(nG, dirla, ed):
    unfound = []
    unfound_rel = []
    fileindex = {}
    fcount = 0
    ccount = 0
    for node in nG.nodes():
        ccount += 1
        if ccount % 20 == 0:
            mfio.save_matrix(fileindex, dirla + 'fileindex_'
                             + str(ccount) + '.dat')
        if nG.node[node]['type'] == 'class':
            if ed + ':-:' + node in fileindex:
                continue
            fileindex[ed + ':-:' + node] = 'inst_' + str(fcount)
            fcount += 1
            instances = ret_inst(ad[ed], node)
            print (node, ed, ccount, len(instances))
            if len(instances) == 0:
                unfound.append((node, ed))
                continue
            f = open(dirla + fileindex[ed + ':-:' + node] + '.tsv', 'w+'
                     )
            f.write(ed + '\t' + node + '\n')
            for k in instances:
                try:
                    f.write(k.encode('utf-8') + '\n')
                except:
                    continue
            f.close()
        else:
            for c in nG.in_edges(node):
                try:
                    if ed + ':-:' + node + ':-:' + c[0] in fileindex:
                        continue
                    fileindex[ed + ':-:' + node + ':-:' + c[0]] = \
                        'rel_' + str(fcount)
                    fcount += 1
                    relations = ret_rels(ad[ed], c[0], c[1])
                    print (node, c[0], ed, ccount, len(relations))
                    if len(relations) == 0:
                        unfound_rel.append((node, c[0], ed))
                        continue
                    f = open(dirla + fileindex[ed + ':-:' + node + ':-:'
                              + c[0]] + '.tsv', 'w+')
                    f.write(ed + '\t' + node + '\t' + c[0] + '\n')
                    for k in relations:
                        try:
                            f.write(k.encode('utf-8') + '\n')
                        except:
                            continue
                    f.close()
                except:
                    continue
    mfio.save_matrix(fileindex, dirla + 'fileindex.dat')
    mfio.save_matrix(unfound, dirla + 'unfound.dat')
    mfio.save_matrix(unfound_rel, dirla + 'unfound_rel.dat')
    return (fileindex, unfound, unfound_rel)


ad = {}
ad["EBI"] = "https://www.ebi.ac.uk/rdf/services/sparql"
ad["PathwayCommons"] = "http://rdf.pathwaycommons.org/sparql/"
ad["LinkedSPL"] = "http://dbmi-icode-01.dbmi.pitt.edu/sparql"
ad["DisGenet"] = "http://rdf.disgenet.org/sparql/"
ad["NextProt"] = "https://api.nextprot.org/sparql"
ad["WikiPathways"] = "http://sparql.wikipathways.org/"
ad["ModelOrganism"] = "http://vt.mo-ld.org/sparql"
ad["EBI-Uniprot"] = "http://sparql.uniprot.org/sparql"
ad["MannheimRDF"] = "http://localhost:8890/sparql"
ad["Bio2Rdf"] = "http://sparql.openlifedata.org/"
ad["NBDC"] = "http://integbio.jp/rdf/sparql"
ad["SBG-TCGA"] = "https://opensparql.sbgenomics.com/blazegraph/namespace/tcga_metadata_kb/sparql"
ad["SBG-CCLE"] = "https://opensparql.sbgenomics.com/blazegraph/namespace/ccle_metadata_kb/sparql"
ad["Mesh"] = "http://id.nlm.nih.gov/mesh/sparql"
ad["OHDSI"] = "http://virtuoso.ohdsi.org:8890/sparql"
ad["Pubchem"] = "https://pubchemdocs.ncbi.nlm.nih.gov/rdf"
ad["linkeddrugs"] = "http://linkeddata.finki.ukim.mk/sparql"
ad["Linkedlifedata"] = "http://linkedlifedata.com/sparql"
ad["LinkedTCGA"] = "http://tcga.deri.ie/"
fdr = 'inst_exp3/'

comG = nx.read_gpickle('lod_schemas3/combG_src.gpickle')
mfio = MatrixIO()

for k in ad:
    nodeset = []
    if k == 'Bio2RDF':
        for node in comG.nodes():
            if 'dsets' not in comG.node[node]:
                continue
            for source in comG.node[node]['dsets'].split(':-:'):
                if source not in ad:
                    nodeset.append(node)
    else:
        for node in comG.nodes():
            if 'dsets' not in comG.node[node]:
                continue
            if k in comG.node[node]['dsets'].split(':-:'):
                nodeset.append(node)
    nG = nx.subgraph(comG, nodeset)
    print 'Now processing ' + k + ' with nodeset ' + str(len(nodeset)) \
        + ' and edge set ' + str(len(nG.edges()))
    dirla = fdr + k + '/'
    subprocess.call(['mkdir', dirla])
    (fileindex, unfound, unfound_rel) = query_instances(nG, dirla, k)
    print 'Processed ' + k + ' with nodeset ' + str(len(fileindex)) \
        + ' and unfound (' + str(len(unfound)) + ', ' \
        + str(len(unfound_rel)) + ')'
