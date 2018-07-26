import json, re, urllib2, sys, os
#import hydra.tpf
import rdflib, hashlib
import cPickle as pickle
import numpy as np
import pandas as pd
from copy import deepcopy

class SPARQLOperator(object):
	def __init__(self, config_filename=None, query_filename=None):
		if config_filename:
			self.config_details = self.read_json(config_filename)
			self.prefix_header = self.set_prefix_header(self.config_details["prefixes"])
		if query_filename:
			self.queries = self.read_json(query_filename)

	def read_json(self, filename):
		with open(filename) as config:
			details = json.load(config)
		return details

	def getPredicate(self, term):
		predicateParts = re.split('[^0-9a-zA-Z_\-]+', term)
		return predicateParts[len(predicateParts)-1]

	def getId(self, uri):
		uriParts = uri.split(":")
		return uriParts[len(uriParts)-1]

	def getLabel(self, label):
		labelParts = label.split("@")
		if len(labelParts) > 1:
			ref_label = labelParts[0][1:-1]
			return ref_label
		else:
			return label

	def set_prefix_header(self, prefixes):
		prefix_header = ""
		for prefix in prefixes:
			prefix_header = prefix_header + "PREFIX " + prefix + ": <" + prefixes[prefix] + "> "
		return prefix_header

	def get_endpoint(self, datasource, endpoints):
		# maybe modify this
		g = rdflib.Graph('TPFStore')
		g.open(endpoints[datasource]["link"])
		return g

	def get_simple_results(self, ontology, query_type):
		g = self.get_endpoint(ontology, self.config_details["endpoints"])
		query =  self.prefix_header + "".join(self.queries[ontology][query_type])
		results = g.query(query)
		return results

# alternatives to creating an object?
class MatrixIO(object):
	def save_matrix(self, x, filename):
		with open(filename, 'wb') as outfile:
			pickle.dump(x, outfile, pickle.HIGHEST_PROTOCOL)

	def load_matrix(self, filename):
		with open(filename, 'rb') as infile:
			x = pickle.load(infile)
		return x

	def save_matrix_txt(self, x, filename):
		np.savetxt(filename, x, delimiter="\t")

	def load_matrix_txt(self, filename):
		return np.matrix(np.loadtxt(filename, delimiter="\t"))

class FileUtils(object):
	def get_reqd_fileset(self, folder_path, exclusion_criteria=None):
		"A class of navigators on the file system"
		fileset = None
		for root, dirs, files in os.walk(folder_path):
			fileset = [x for x in files if not exclusion_criteria(x)] if exclusion_criteria else files
		return fileset

	def assign_filesize(self, folder_path, exclusion_criteria=None):
		'''Returns a dict with file name as key and its size (bytes) as value'''
		fileset = self.get_reqd_fileset(folder_path, exclusion_criteria)
		file_dict = {}
		for _file in fileset:
			size = os.stat(folder_path + _file).st_size
			file_dict[_file] = size
		return file_dict
		'''returns a simple dict indexed by values in the given column from a TSV file'''
	'''
	def get_simple_dict(self, file_name, column=0, header=None, names=[], cols=[]):
		f = open(file_name)
		flines = f.readlines()
		f.close()
		fdict = {}
		for k in range(len(flines)):
			if header and k==0:
				continue
			fparts = flines[k].strip().split("\t")
			fdict[fparts[column]] = {names[m]: fparts[cols[m]] for m in cols}
		return fdict'''


class DictUtils(object):
	def merge_two_dicts(self, x, y):
		'''Given two dicts, merge them into a new dict as a shallow copy.'''
		z = x.copy()
		z.update(y)
		return z

	def iterative_merge_pickles(self, folder_path, filename_generator=None, exclusion_criteria=None, final_name=None, verbose=True, total_iter=None):
		'''Iteratively merge dictionaries stored as pickles in the folder_path, adhering to certain exclusion_criteria (excrit=True if you want to exclude those files'''
		mfio = MatrixIO()
		found_true_love = False
		num_iter = 0
		total_iter = total_iter if total_iter else 10
		while True:
			# Use a queue data structure here ... one pop for true
			if found_true_love or num_iter > total_iter:
				# rename file to final_name
				break
			num_iter += 1
			fu = FileUtils()
			fileset = fu.get_reqd_fileset(folder_path, exclusion_criteria)
			print (fileset)
			exists_current = None
			if len(fileset) == 1:
				found_true_love = True
				break
			for file in fileset:
				if not exists_current:
					exists_current = file
				else:
					dict1 = mfio.load_matrix(folder_path + "/" + file)
					dict2 = mfio.load_matrix(folder_path + "/" + exists_current)
					merged_dict = self.merge_two_dicts(dict1, dict2)
					if not filename_generator:
						# this is a last minute resort assuming the user is too lazy to provide a filename_generator function
						filename = hashlib.sha1(exists_current + "-" + file) 
					else:
						filename = filename_generator(exists_current, file)
					mfio.save_matrix(merged_dict, folder_path + "/" + filename)
					os.remove(folder_path + "/" + exists_current)
					os.remove(folder_path + "/" + file)
					if verbose:
						print ("Merged file " + exists_current + " and file " + file + " into " + filename)
					exists_current = None 

class FrameUtils(object):
	def convertdf(self, df, column):
		''' this is a new way to convert a dataframe into a dictionary, such that values of a column end up becoming the keys in the dictionary
		This function should be moved to Utils'''
		# such that values in one column become the keys
		dfc = deepcopy(df)
		dfc = dfc.set_index(column) 
		return dfc.to_dict(orient="index")


class HTTPUtils(object):
	def get_json(self, url, headers=None):
	    try:
	        opener = urllib2.build_opener()
	        if headers:
		        opener.addheaders = headers
	        return json.loads(opener.open(url).read())
	    except urllib2.HTTPError:
	        print ('Error Encountered')
	        return None

def test_dictutils():
	du = DictUtils()
	def expcri(x):
		# this should be your own exclusion criteria
		#print x
		xparts = x.split(".")
		if xparts[1] != "dat" or xparts[0][0:5] != "assoc":
			return True
		return False
	def filename_gen(name1, name2):
		name1parts = name1.split(".")[0].split("_")
		name2parts = name2.split(".")[0].split("_")
		new_name = "_".join(name1parts[0:len(name1parts)-1]) + "_"
		new_name = new_name + hashlib.sha1(name1parts[len(name1parts)-1] + "-" + name2parts[len(name2parts)-1]).hexdigest()
		new_name = new_name + ".dat"
		return new_name
	du.iterative_merge_pickles("dat_files", filename_generator=filename_gen, exclusion_criteria=expcri, final_name="assoc_paths_final.dat", verbose=True)

def test_httputils():
	hu = HTTPUtils()
	REST_URL = "http://data.bioontology.org"
	API_KEY = "5b7f7e20-015c-496f-ade9-ca341345cef7"
	headers = [('Authorization', 'apikey token=' + API_KEY)]
	# = hu.get_json(REST_URL+"/annotator", headers)

def merge(x):
	"merge lists of list in one list"
	a = []
	for m in x:
		a.extend(m)
	return a

if __name__ == "__main__":
	test_dictutils()

	
