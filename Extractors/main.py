import json, re, sys, os, time
from operator import itemgetter
import numpy as np
import pandas as pd
import networkx as nx
from copy import deepcopy
from utils import HTTPUtils, MatrixIO, FileUtils

first_cap_re = re.compile('(.)([A-Z][a-z]+)')
all_cap_re = re.compile('([a-z0-9])([A-Z])')

def parse_camel_case(name):
    s1 = first_cap_re.sub(r'\1_\2', name)
    return all_cap_re.sub(r'\1_\2', s1).lower()

def parse_uri_token(token):
    prop = parse_camel_case(token)
    pparts = re.split("[-_]", prop)
    name = " ".join([x.title() for x in pparts])
    return name

def parseURI(uri):
    nparts = re.split("[#:/]", uri)
    name = parse_uri_token(nparts[len(nparts)-1])
    return name

def model_proc(k):
    if "resource" in k:
        a = parseURI(k).split()
        if len(a) > 1: return "ModelOrganism-" + a[0]
    return "ModelOrganism"

def nbdc_proc(node):
    nparts = re.split("[#:/]", node)
    prefix = node[:len(node)-len(nparts[len(nparts)-1])]
    if prefix in ndbc: return ndbc[prefix]["id"]
    else: return "NBDC"

def bio2rdf_proc(k):
    dsets = []
    if k in ap:
        for m in ap[k]:
            if "bio2rdf" in m:
                a = m.split(".")
                dset = "Bio2Rdf-" + a[len(a)-2].title()
                dsets.append(dset)
            elif "disgenet" in m:
                dset = "DisGenet"
                dsets.append(dset)
        return ":-:".join(dsets)
    else: return "Bio2Rdf"

def pref_dsets_cst(dset_l):
    ndsets = []
    dsets = dset_l.split(":-:")
    for m in dsets:
        if m == 'Pdb': ndsets.append("NBDC")
        elif m == 'BioSamples': ndsets.append('EBI-BioSamples')
        elif m == 'tcga': ndsets.append("LinkedTCGA")
        else: ndsets.append(m)
    return ":-:".join(ndsets)

def det_dset(ami):
    (a, m, i) = ami
    dsets = []
    inst_count = []
    for xi in range(len(m)): 
        x = m[xi]
        ic = i[xi]
        nist = False
        if "linkedlifedata" in x.lower(): 
            dsets.append(x.split("_")[2].split(".")[0].title())
        elif "pubchem" in x.lower(): 
            dsets.append("Pubchem")
        elif "modelorganism" in x.lower(): 
            dsets.append(model_proc(a))
        elif "bio2rdf-" in x.lower(): 
            dsets.append(x.split("_")[2].split(".")[0].title())
        elif "bio2rdf" in x.lower(): 
            b2r_dsets = bio2rdf_proc(a)
            dsets.append(b2r_dsets)
            nist = True
            for k in range(len(b2r_dsets.split(":-:"))):
                inst_count.append(ic)
        elif "nbdc" in x.lower(): 
            dsets.append(nbdc_proc(a))
        else: 
            dsets.append(x.split("_")[2].split(".")[0])
        if not nist:
            inst_count.append(ic)
    dsets = ":-:".join([str(k) for k in dsets])
    return dsets, inst_count

def get_sparql_endpoint(dset):
    if dset in ad: 
        return ad[dset]
    elif dset.split("-")[0] in ad: 
        return ad[dset.split("-")[0]]
    else: 
        return ""

def get_folder_path(dset):
    folder_path = final_folder + "instance_data/"
    if "bio2rdf" in dset.lower():
        subdir = os.listdir(folder_path + "Bio2Rdf/")
        if dset in subdir: return folder_path + "Bio2Rdf/" + dset + "/"
        else: return folder_path + "Bio2Rdf/Bio2Rdf/"
    elif "EBI" in dset:
        return folder_path + "EBI/" + dset + "/"
    elif "sbg" in dset.lower():
        return folder_path + "SBG/" + dset + "/"
    elif "nbdc" in dset.lower():
        return folder_path + "NBDC/"
    elif "linkedlifedata" in dset.lower():
        return folder_path + "Linkedlifedata/" + dset + "/"
    elif "modelorganism" in dset.lower():
        return folder_path + "ModelOrganism/"
    else:
        return folder_path + dset + "/"
    
def get_inst_data(dset, verbose=True):
    str_inst_val = True
    num_inst_val = True
    if len(dset) == 0: 
        return None, None, False, False
    folder_path = get_folder_path(dset)
    if verbose: print ("---------------------------------")
    if verbose: print (dset, folder_path)
    if os.path.isdir(folder_path):
        str_inst = pd.read_csv(folder_path + "str_inst_df.tsv", sep="\t")
        num_inst = pd.read_csv(folder_path + "num_inst_df.tsv", sep="\t")
        if verbose: print ("Dset sizes", str_inst.shape, num_inst.shape)
        if str_inst.shape[0] == 0: str_inst_val = False
        if num_inst.shape[0] == 0: num_inst_val = False
        return str_inst, num_inst, str_inst_val, num_inst_val
    else: 
        print ("No folder path found -  XXX")
        return None, None, None, None

def get_new_dsets(node):
    src_dsets = combG.node[node]["dsets"] if "dsets" in combG.node[node] else ""
    src_dsets = [old_dset_mappings[k] for k in src_dsets.split(":-:") if k in old_dset_mappings]
    return src_dsets

def get_assoc_properties(source, dset):
    data_properties = []
    obj_properties = []
    def_l_set = 'http://www.w3.org/2001/XMLSchema#string'
    if not combG.has_node(source): return ([], [])
    src_dsets = get_new_dsets(source)
    #print (src_dsets)
    if "NBDC" in dset: dset = "NBDC"
    elif "ModelOrganism" in dset: dset = "ModelOrganism"
    else: dset = dset
    for k in combG[source]:
        dprop_material = set([])
        oprop_material = set([])
        if combG.node[k]["type"] == "literal" or combG.node[k]["type"] == "objectProp":
            targ_dsets = get_new_dsets(k)
            dset_intersection = set(src_dsets).intersection(set(targ_dsets))
            if dset in dset_intersection:
                for n in combG[source][k]:
                    mat = combG[source][k][n]
                    _type = mat['type'] if 'type' in mat else def_l_set
                    _count = mat['count'] if 'count' in mat else 0
                    if len(_type.split("->")) == 1: 
                        dprop_material.add((_count, _type))
                    else:
                        oprop_material.add((_count, _type))
                if len(dprop_material) > 0: 
                    data_properties.append((k, dprop_material))
                if len(oprop_material) > 0:
                    obj_properties.append((k, oprop_material))
    return (data_properties, obj_properties)

def get_labels(node):
    parsed_label = combG.node[node]["refL1"] if "refL1" in combG.node[node] else parseURI(node)
    sec_label = combG.node[node]["sec_label"] if "sec_label" in combG.node[node] else ""
    return parsed_label, sec_label

def prepare_property_dtype_dict(dset):
    curop = {}
    curdp = {}
    dtypes = set([])
    opn_c = 0
    dpn_c = 0
    for m in all_dsets[dset]['Classes']:
        dprop, oprop = get_assoc_properties(m, dset)
        for item in oprop:
            n = item[0]
            mn = item[1]
            if not n in curop: curop[n] = {"Realizations": {}}
            parsed_label, sec_label = get_labels(n)
            curop[n]["ParsedLabel"] = parsed_label
            curop[n]["ExternalLabel"] = sec_label
            curop[n]["Shared Datasets"] = get_new_dsets(n)
            if n in oproponto_nodes: curop[n]["FromBioOntology"] = oproponto_nodes[n]
            else: curop[n]["FromBioOntology"] = ""
            if n in opropvocab_nodes: curop[n]["FromLOVVocabulary"] = opropvocab_nodes[n]
            else: curop[n]["FromLOVVocabulary"] = ""
            for op in mn:
                opd = op[1].split("->")
                _domain = m
                _range = opd[2] if len(opd) > 2 else ""
                _count = op[0]
                pdict = {"Domain": _domain, "Range": _range, "Assertion Count": _count}
                curop[n]["Realizations"][op[1]] = pdict
                opn_c += 1
        for item in dprop:
            n = item[0]
            mn = item[1]
            if not n in curdp: curdp[n] = {"Realizations": {}}
            parsed_label, sec_label = get_labels(n)
            curdp[n]["ParsedLabel"] = parsed_label
            curdp[n]["ExternalLabel"] = sec_label
            curdp[n]["Shared Datasets"] = get_new_dsets(n)
            if n in dproponto_nodes: curdp[n]["FromBioOntology"] = dproponto_nodes[n]
            else: curdp[n]["FromBioOntology"] = ""
            if n in dpropvocab_nodes: curdp[n]["FromLOVVocabulary"] = dpropvocab_nodes[n]
            else: curdp[n]["FromLOVVocabulary"] = ""
            for dp in mn:
                _range = dp[1]
                _count = dp[0]
                dtypes.add(_range)
                try:
                    realization = m + "->" + n + "->" + _range
                    pdict = {"Domain": m, "Range": _range, "Assertion Count": _count}
                    curdp[n]["Realizations"][realization] = pdict
                    dpn_c += 1
                except: 
                    print ("Error generating realization for", dp, n, m)
    dtypes = list(dtypes)
    return curop, curdp, dtypes, opn_c, dpn_c

def get_property_sets_refine(psets, ponto_nodes, pvocab_nodes, ptype):
    property_sets = deepcopy(psets)
    for k in property_sets:
        parsed_label, sec_label = get_labels(k)
        #print (parsed_label, sec_label)
        if k in ponto_nodes: oin_onto = ponto_nodes[k]
        else: oin_onto = ""
        if k in pvocab_nodes: oin_vocab = pvocab_nodes[k]
        else: oin_vocab = ""
        shared_dsets = []
        realization_count = []
        assertion_count = []
        for m in all_dsets:
            if k in all_dsets[m]["Properties"][ptype]:
                shared_dsets.append(m)
                real = all_dsets[m]["Properties"][ptype][k]['Realizations']
                realization_count.append(int(len(real)))
                assertions = int(np.sum([real[n]['Assertion Count'] for n in real]))
                assertion_count.append(assertions)
        property_sets[k] = {"ParsedLabel": parsed_label, "ExternalLabel": sec_label, "FromBioOntology": oin_onto,
                            "FromLOVVocabulary": oin_vocab, "Datasets": ":-:".join(shared_dsets), 
                            "DatasetCount": len(shared_dsets), 
                            "RealizationCount": ":-:".join([str(n) for n in realization_count]), 
                            "TotalRealizationCount": np.sum(realization_count), 
                            "AssertionCount": ":-:".join([str(n) for n in assertion_count]), 
                            "TotalAssertionCount": np.sum(assertion_count)}
    property_sets_df = pd.DataFrame.from_dict(property_sets, orient="index")
    property_sets_df = property_sets_df.reset_index()
    property_sets_df = property_sets_df.sort_values("TotalAssertionCount", ascending=False)
    property_sets_df = property_sets_df[[u'index', u'ParsedLabel', u'ExternalLabel', u'FromBioOntology', 
                                u'FromLOVVocabulary', u'Datasets', u'DatasetCount', 
                                u'RealizationCount', u'TotalRealizationCount', u'AssertionCount', 
                                u'TotalAssertionCount']]
    return property_sets_df

print ("Reading main data files")
final_folder = "raw_data_folder/"
ap = mfio.load_matrix(final_folder + "network_data/bio2rdf_classdsets.dat")
mfio = MatrixIO()
fu = FileUtils()
combG = nx.read_gpickle(final_folder + "network_data/combG_c123.gpickle")
ref_preflist = mfio.load_matrix(final_folder + "network_data/ref-preflist.dat")
ref_obolist = mfio.load_matrix(final_folder + "network_data/ref-oboprefset.dat")
vocab_elem = mfio.load_matrix(final_folder + "network_data/vocab_elem.dat")
sc_nomen = mfio.load_matrix(final_folder + "network_data/sc_nomen.dat")
print ("Vocabularies from LOV", len(sc_nomen))

print ("Reading auxiliary data files")

ad = pd.read_csv(final_folder + "util_files/LDProjects.tsv", sep="\t").set_index("index")
ad = ad.to_dict(orient="index")
for k in ad: ad[k] = ad[k]['0']
print ("Linked Data Projects", len(ad))

ndbc = pd.read_csv(final_folder + "util_files/NBDC_prefix.tsv", sep="\t").set_index("index")
ndbc = ndbc.to_dict(orient="index")
print ("NBDC Graphs", len(ndbc))

disAllPrefs = pd.read_csv(final_folder + "util_files/ld_onto_prefixes.tsv", sep="\t").set_index("index")
disAllPrefs = disAllPrefs.to_dict(orient="index")
for k in disAllPrefs: disAllPrefs[k] = disAllPrefs[k]['0']
print("Linked data ontology prefixes", len(disAllPrefs))

ontoLabels = pd.read_csv(final_folder + "util_files/ontology_labels.tsv", sep="\t").set_index("index")
ontoLabels = ontoLabels.to_dict(orient="index")
print("Ontology Labels", len(ontoLabels))

old_dset_mappings = pd.read_csv(final_folder + "util_files/old_dataset_mappings.tsv", sep="\t").set_index("index")
old_dset_mappings = old_dset_mappings.to_dict(orient="index")
for k in old_dset_mappings: old_dset_mappings[k] = old_dset_mappings[k]['0']
old_dset_mappings[''] = ''
print ("Old Dataset Mappings", len(old_dset_mappings))

class_list_folder = final_folder + "class_lists/"
class_set_files = fu.get_reqd_fileset(class_list_folder, lambda x: False if "classlist" in x.lower() else True)
class_sets = {}
ccount = 0
for k in class_set_files:
    a = mfio.load_matrix(class_list_folder + k)
    for m in a: 
        if not m in class_sets: class_sets[m] = {"files": [], "instance_count": []}
        class_sets[m]["files"].append(k)
        class_sets[m]["instance_count"].append(a[m])
        ccount += 1

print ("Total Class Count", ccount)
print ("Unique Class Count", len(class_sets))
for k in class_sets:
    a, b = det_dset((k, class_sets[k]["files"], class_sets[k]["instance_count"]))
    class_sets[k]["dsets"] = a
    class_sets[k]["instance_count"] = b

for k in class_sets:
    class_sets[k]["total_instance_count"] = np.sum(class_sets[k]["instance_count"])

uc = []
undsets = {}
for k in combG.nodes():
    if not combG.node[k]["type"] == "class": continue
    if not k in class_sets:
        uc.append((k, combG.node[k]))
        if "dsets" in combG.node[k]: 
            spcdsets = combG.node[k]["dsets"].split(":-:")
            for m in spcdsets:
                if not m in undsets: undsets[m] = []
                undsets[m].append(k)

for k in undsets:
    for m in undsets[k]:
        class_sets[m] = {"files": [], 
                         "dsets": pref_dsets_cst(combG.node[m]["dsets"]), 
                         "instance_count": [combG.node[m]["count"]] if "count" in combG.node[m] else [0], 
                         "total_instance_count": combG.node[m]["count"] if "count" in combG.node[m] else 0}

found_lab = 0
found_sec_lab = 0
for k in class_sets:
    if combG.has_node(k):
        if "refL1" in combG.node[k]: 
            class_sets[k]["refL1"] = combG.node[k]["refL1"]
            found_lab += 1
        if "sec_label" in combG.node[k]:
            class_sets[k]["sec_label"] = combG.node[k]["sec_label"]
            found_sec_lab += 1
    if not "refL1" in class_sets[k]: class_sets[k]["refL1"] = parseURI(k)
    if not "sec_label" in class_sets[k]: class_sets[k]["sec_label"] = ""

print ("Found parsed labels", found_lab)
print ("Found external labels", found_sec_lab)

sort_obolist = sorted({k: len(ref_obolist[k]) for k in ref_obolist}.items(), key=itemgetter(1), reverse=True)
for k in sort_obolist: disAllPrefs[k[0]] = parseURI(k[0]).upper()

onto_nodes = {}
for k in disAllPrefs:
    if k in ref_obolist: ic_n = ref_obolist[k]
    else: ic_n = ref_preflist[k]
    for m in ic_n:
        if m not in class_sets: continue
        onto_nodes[m] = ontoLabels[k]['name']
print ("Ontology related nodes (only classes)", len(onto_nodes))

vocab_index = {}
for k in vocab_elem:
    for m in vocab_elem[k]:
        if m not in vocab_index: vocab_index[m] = set([])
        vocab_index[m].add(k)
print ("Terms from vocabularies", len(vocab_index))

vocab_nodes = {}
for k in vocab_index:
    st = list(vocab_index[k])[0]
    vocabiden = sc_nomen[st]["vocab_name"].title()
    if not k in onto_nodes:
        if k in class_sets: vocab_nodes[k] = vocabiden
print ("Vocabulary related nodes (only classes)", len(vocab_nodes))

for k in class_sets:
    if k in onto_nodes: class_sets[k]["FromOntology"] = onto_nodes[k]
    else: class_sets[k]["FromOntology"] = ""
    if k in vocab_nodes: class_sets[k]["FromVocabulary"] = vocab_nodes[k]
    else: class_sets[k]["FromVocabulary"] = ""

uniq_classes = []
for k in class_sets:
    if class_sets[k]["FromOntology"] == "" and class_sets[k]["FromVocabulary"] == "":
        uniq_classes.append(k)
        
print ("Unique non-reused classes", len(uniq_classes))

class_sets_df = pd.DataFrame.from_dict(class_sets, orient="index")
class_sets_df = class_sets_df.reset_index()
del class_sets_df["files"]
class_sets_df = pd.concat([class_sets_df[["index", "refL1", "sec_label", "FromOntology", "FromVocabulary", "dsets"]], 
                           class_sets_df["dsets"].apply(lambda x: len(x.split(":-:"))).to_frame(name="dcount"), 
                           class_sets_df["instance_count"].apply(lambda x: ":-:".join([str(k) for k in x])), 
                           class_sets_df["total_instance_count"]], axis=1)
class_sets_df = class_sets_df.sort_values(["total_instance_count"], ascending=False)
class_sets_df.columns = ["URI", "ParsedLabel", "ExternalLabel", "FromBioOntology", "FromLOVVocabulary", "Datasets", 
                         "DatasetCount", "InstanceCountByDataset", "TotalInstanceCount"]

print ("Total Classes", class_sets_df.shape[0])

class_sets_df.to_csv("final/Extracted_Classes.tsv", sep="\t", index=None)

dproperty_sets = {}
oproperty_sets = {}
for k in combG.nodes():
    if not "type" in combG.node[k]: continue
    if combG.node[k]["type"] == 'literal': dproperty_sets[k] = {}
    if combG.node[k]["type"] == 'objectProp': oproperty_sets[k] = {}

print ("Object Properties", len(oproperty_sets))
print ("Data Properties", len(dproperty_sets))

oproponto_nodes = {}
dproponto_nodes = {}
for k in disAllPrefs:
    if k in ref_obolist: ic_n = ref_obolist[k]
    else: ic_n = ref_preflist[k]
    for m in ic_n:
        if m in oproperty_sets:
            oproponto_nodes[m] = ontoLabels[k]['name']
        elif m in dproperty_sets:
            dproponto_nodes[m] = ontoLabels[k]['name']
print ("Ontology related nodes (only object properties)", len(oproponto_nodes))
print ("Ontology related nodes (only data properties)", len(dproponto_nodes))

opropvocab_nodes = {}
dpropvocab_nodes = {}
for k in vocab_index:
    st = list(vocab_index[k])[0]
    vocabiden = sc_nomen[st]["vocab_name"].title()
    if not k in oproponto_nodes:
        if k in oproperty_sets: opropvocab_nodes[k] = vocabiden
    if not k in dproponto_nodes:
        if k in dproperty_sets: dpropvocab_nodes[k] = vocabiden
print ("Vocabulary related nodes (only object properties)", len(opropvocab_nodes))
print ("Vocabulary related nodes (only data properties)", len(dpropvocab_nodes))

all_dsets = {}
for k in set(class_sets_df["Datasets"]):
    ap = k.split(":-:")
    for m in ap: 
        if not m in all_dsets: 
            all_dsets[m] = {"SPARQL Endpoint": get_sparql_endpoint(m), 
                            "Counts": {
                                "Class Count": 0,
                                "Object Property Count": 0,
                                "Data Property Count": 0,
                                "Datatypes Count": 0,
                                "Object Property Realization Count": 0,
                                "Data Property Realization Count": 0,
                                "Class-with-Instance Count": 0,
                                "Object Property-with-Assertion Count": 0,
                                "Data Property-with-Assertion Count": 0
                            },
                            "Classes": {},
                            "Properties": {
                                "Object Properties": {},
                                "Data Properties": {}
                            },
                            "Datatypes": []
                           }
    
print ("Total Datasets", len(all_dsets))

for k in class_sets:
    dsets_c = class_sets[k]['dsets'].split(":-:")
    instances = class_sets[k]['instance_count']
    for m in range(len(dsets_c)):
        class_dict = {"URI": k, 
                      "Parsed Label": class_sets[k]["refL1"], 
                      "External Label": class_sets[k]["sec_label"], 
                      "Shared Datasets": class_sets[k]["dsets"].split(":-:"), 
                      "Source BioPortal Ontology": class_sets[k]["FromOntology"],
                      "Source LOV Vocabulary": class_sets[k]["FromVocabulary"],
                      "Instance Count": instances[m], 
                      "Instance Characteristics": {}}
        all_dsets[dsets_c[m]]['Classes'][k] = class_dict
        all_dsets[dsets_c[m]]["Counts"]['Class Count'] += 1

for k in all_dsets:
    str_inst, num_inst, str_valid, num_valid = get_inst_data(k)
    data_count = 0
    if str_valid:
        entity_data = str_inst[str_inst['ftype'] == 'entity']
        entity_data = entity_data.set_index("index").to_dict(orient="index")
        for m in entity_data:
            if m in all_dsets[k]['Classes']:
                instance_dict = {"Sample Instances (3)": eval(entity_data[m]['samp_insts']), 
                                 "Namespaces": list(eval(entity_data[m]['namespace'])), 
                                 "Is Categorical": entity_data[m]['is_categorical'], 
                                 "Sample Instances Total": entity_data[m]['inst_count'], 
                                 'Categories': list(eval(entity_data[m]['categories'])), 
                                 'Datatype': entity_data[m]['dtype'], 
                                 'Instance Length (Median)': entity_data[m]['lenmedian'], 
                                 'Instance Length (Standard Deviation)': entity_data[m]['lenstd'], 
                                 'Is String': True}
                all_dsets[k]['Classes'][m]["Instance Characteristics"] = instance_dict
                all_dsets[k]["Counts"]["Class-with-Instance Count"] += 1

for k in sorted(all_dsets.keys()):
    odict, ddict, datatypes, opn, dpn = prepare_property_dtype_dict(k)
    print (k, len(odict), len(ddict), len(datatypes))
    all_dsets[k]["Properties"]["Object Properties"] = odict
    all_dsets[k]["Properties"]["Data Properties"] = ddict
    all_dsets[k]["Datatypes"] = datatypes
    all_dsets[k]["Counts"]["Object Property Count"] = len(odict)
    all_dsets[k]["Counts"]["Data Property Count"] = len(ddict)
    all_dsets[k]["Counts"]["Datatypes Count"] = len(datatypes)
    all_dsets[k]["Counts"]["Object Property Realization Count"] = opn
    all_dsets[k]["Counts"]["Data Property Realization Count"] = dpn

found_op_s = 0
found_op_n = 0
found_dp_s = 0
found_dp_n = 0

for k in all_dsets:
    op = all_dsets[k]["Properties"]['Object Properties']
    dp = all_dsets[k]["Properties"]['Data Properties']
    str_inst, num_inst, str_valid, num_valid = get_inst_data(k)
    data_count = 0
    str_relation_data = {}
    num_relation_data = {}
    if str_valid:
        str_relation_data = str_inst[str_inst['ftype'] == 'relation']
        str_relation_data = str_relation_data.set_index("index").to_dict(orient="index")
    if num_valid:
        num_relation_data = num_inst[num_inst['ftype'] == 'relation']
        num_relation_data = num_relation_data.set_index("index").to_dict(orient="index")
    for rel in op:
        has_realization = False
        for dreal in op[rel]['Realizations']:
            dreal_sp = "->".join(dreal.split("->")[0:2])
            if dreal_sp in str_relation_data:
                try:
                    instance_dict = {"Sample Instances (3)": eval(str_relation_data[dreal_sp]['samp_insts']), 
                                     "Namespaces": list(eval(str_relation_data[dreal_sp]['namespace'])), 
                                     "Is Categorical": str_relation_data[dreal_sp]['is_categorical'], 
                                     "Sample Instances Total": str_relation_data[dreal_sp]['inst_count'],
                                     'Categories': list(eval(str_relation_data[dreal_sp]['categories'])), 
                                     'Datatype': str_relation_data[dreal_sp]['dtype'], 
                                     'Instance Length (Median)': str_relation_data[dreal_sp]['lenmedian'],
                                     'Instance Length (Standard Deviation)': str_relation_data[dreal_sp]['lenstd'], 
                                     'Is String': True}
                    op[rel]['Realizations'][dreal]["Assertion Characteristics"] = instance_dict
                    found_op_s += 1
                except: continue
            elif dreal_sp in num_relation_data:
                try:
                    instance_dict = {"Sample Instances (3)": eval(num_relation_data[dreal_sp]['samp_insts']), 
                                     "Is Categorical": num_relation_data[dreal_sp]['is_categorical'], 
                                     "Sample Instances Total": num_relation_data[dreal_sp]['inst_count'],
                                     'Categories': list(eval(num_relation_data[dreal_sp]['categories'])), 
                                     'Datatype': num_relation_data[dreal_sp]['dtype'], 
                                     'Value Set (Median)': num_relation_data[dreal_sp]['valmedian'],
                                     'Value Set (Standard Deviation)': num_relation_data[dreal_sp]['valstd'], 
                                     'Normal Distribution Test (P-Value)': num_relation_data[dreal_sp]['normp'],
                                     'Normal Distribution Test (Statistic)': num_relation_data[dreal_sp]['normsk'],
                                     'Is String': False}
                    op[rel]['Realizations'][dreal]["Assertion Characteristics"] = instance_dict
                    found_op_n += 1
                except: continue
            op[rel]['Realizations'][dreal]['Assertion Count'] = int(op[rel]['Realizations'][dreal]['Assertion Count'])
            if op[rel]['Realizations'][dreal]['Assertion Count'] > 0: has_realization = True
        if has_realization:
            all_dsets[k]["Counts"]["Object Property-with-Assertion Count"] += 1
    for rel in dp:
        has_realization = False
        for dreal in dp[rel]['Realizations']:
            dreal_sp = "->".join(dreal.split("->")[0:2])
            if dreal_sp in str_relation_data:
                try:
                    instance_dict = {"Sample Instances (3)": eval(str_relation_data[dreal_sp]['samp_insts']), 
                                     "Namespaces": list(eval(str_relation_data[dreal_sp]['namespace'])), 
                                     "Is Categorical": str_relation_data[dreal_sp]['is_categorical'], 
                                     "Sample Instances Total": str_relation_data[dreal_sp]['inst_count'],
                                     'Categories': list(eval(str_relation_data[dreal_sp]['categories'])), 
                                     'Datatype': str_relation_data[dreal_sp]['dtype'], 
                                     'Instance Length (Median)': str_relation_data[dreal_sp]['lenmedian'],
                                     'Instance Length (Standard Deviation)': str_relation_data[dreal_sp]['lenstd'], 
                                     'Is String': True}
                    dp[rel]['Realizations'][dreal]["Assertion Characteristics"] = instance_dict
                    found_dp_s += 1
                except: continue
            elif dreal_sp in num_relation_data:
                try:
                    instance_dict = {"Sample Instances (3)": eval(num_relation_data[dreal_sp]['samp_insts']), 
                                     "Is Categorical": num_relation_data[dreal_sp]['is_categorical'], 
                                     "Sample Instances Total": num_relation_data[dreal_sp]['inst_count'],
                                     'Categories': list(eval(num_relation_data[dreal_sp]['categories'])), 
                                     'Datatype': num_relation_data[dreal_sp]['dtype'], 
                                     'Value Set (Median)': num_relation_data[dreal_sp]['valmedian'],
                                     'Value Set (Standard Deviation)': num_relation_data[dreal_sp]['valstd'], 
                                     'Normal Distribution Test (P-Value)': num_relation_data[dreal_sp]['normp'],
                                     'Normal Distribution Test (Statistic)': num_relation_data[dreal_sp]['normsk'],
                                     'Is String': False}
                    dp[rel]['Realizations'][dreal]["Assertion Characteristics"] = instance_dict
                    found_dp_n += 1
                except: continue
            dp[rel]['Realizations'][dreal]['Assertion Count'] = int(dp[rel]['Realizations'][dreal]['Assertion Count'])
            if dp[rel]['Realizations'][dreal]['Assertion Count'] > 0: has_realization = True
        if has_realization:
            all_dsets[k]["Counts"]["Data Property-with-Assertion Count"] += 1

oproperty_sets_df = get_property_sets_refine(oproperty_sets, oproponto_nodes, opropvocab_nodes, "Object Properties")
dproperty_sets_df = get_property_sets_refine(dproperty_sets, dproponto_nodes, dpropvocab_nodes, "Data Properties")
dproperty_sets_df.to_csv("final/Extracted_Data_Properties.tsv", sep="\t", index=None)
oproperty_sets_df.to_csv("final/Extracted_Object_Properties.tsv", sep="\t", index=None)

all_datatypes = {}
for k in all_dsets:
    dtypes = all_dsets[k]["Datatypes"]
    for m in dtypes:
        if len(m) == 0: continue
        if not m in all_datatypes: all_datatypes[m] = {"Datasets": set([])}
        all_datatypes[m]["Datasets"].add(k)

print (len(all_datatypes))

all_datatypes_df = pd.DataFrame.from_dict(all_datatypes, orient="index")
all_datatypes_df = all_datatypes_df.reset_index()
all_datatypes_df = pd.concat([all_datatypes_df["index"], all_datatypes_df["Datasets"].apply(lambda x: ":-:".join(list(x))), 
                              all_datatypes_df["Datasets"].apply(lambda x: len(x)).to_frame(name="DatasetCount")], axis=1)
all_datatypes_df = all_datatypes_df.sort_values("DatasetCount", ascending=False)
all_datatypes_df.to_csv("final/Extracted_Datatypes.tsv", sep="\t", index=None)
all_datatypes_df.head()

del all_dsets[""]
print ("Removing non-source", len(all_dsets))

dset_print = {}
print_attr = ["SPARQL Endpoint", "Class Count", "Object Property Count", 
              "Data Property Count", "Object Property Realization Count", "Data Property Realization Count", 
              "Datatypes Count", "Class-with-Instance Count", "Data Property-with-Assertion Count", 
              "Object Property-with-Assertion Count"]
for k in all_dsets:
    if len(k) == 0: continue
    dset_print[k] = {}
    for m in print_attr:
        if not "count" in m.lower(): dset_print[k][m] = all_dsets[k][m]
        else: dset_print[k][m] = all_dsets[k]["Counts"][m]
    dset_print[k]["Class Instance Coverage"] = float(all_dsets[k]["Counts"]["Class-with-Instance Count"])/all_dsets[k]["Counts"]["Class Count"]

dset_print = pd.DataFrame.from_dict(dset_print, orient="index")
dset_print = dset_print.reset_index()
dset_print = dset_print[["index", "SPARQL Endpoint", "Class Count", "Object Property Count", 
              "Data Property Count", "Object Property Realization Count", "Data Property Realization Count", 
              "Datatypes Count", "Class-with-Instance Count", 
              "Object Property-with-Assertion Count", "Data Property-with-Assertion Count", "Class Instance Coverage"]]
dset_print.to_csv("final/Linked_Dataset_Graphs.tsv", sep="\t", index=None)

mfio.save_matrix(all_dsets, "final/LSLOD_Schema_Graph.json.pickle")
