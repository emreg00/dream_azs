
from toolbox import configuration, guild_utilities, wrappers
from time import strftime
import os, csv, numpy, cPickle
try:
    from indigo import indigo
except:
    print "Indigo not found, chemical similarity will not be available!" 

CONFIG = configuration.Configuration() 


def main():
    # GUILD based synergy
    #get_guild_based_synergy()
    # Create feature file
    create_feature_file()
    return


def create_feature_file():
    # Get drug info
    drug_to_values = get_drug_info(nodes=None)
    print drug_to_values.items()[:3]
    # Get cell line info
    cell_line_to_value = get_cell_line_info()
    print cell_line_to_value.items()[:3]
    # Get synergy info
    combination_to_values = get_synergy_info()
    print combination_to_values.items()[:3]
    # Get gexp info
    gexp_norm, gene_to_idx, cell_line_to_idx = wrappers.get_expression_info(gexp_file = CONFIG.get("gexp_file"), process=set(["z", "abs"]), dump_file = CONFIG.get("gexp_dump"))
    values = [gene_to_idx["TSPAN6"], gene_to_idx["TNMD"]] 
    print gexp_norm[values, cell_line_to_idx["647-V"]]
    print "TNMD @ 647-V", gexp_norm[gene_to_idx["TNMD"], cell_line_to_idx["647-V"]] 
    # Get mutation info
    gene_to_cell_line_to_mutation = get_mutation_info() 
    print gene_to_cell_line_to_mutation.items()[:3]
    # Get CNV info
    #gene_to_cell_line_to_cnv = {}
    gene_to_cell_line_to_cnv = get_cnv_info(CONFIG.get("cnv_file")) 
    print gene_to_cell_line_to_cnv.items()[:3]
    # Get GUILD info
    #combination_to_guild_values = {}
    combination_to_guild_values = get_guild_based_synergy_scores(drug_to_values.keys(), None, None, None)
    print combination_to_guild_values.items()[:3]
    # Get cancer gene & pathway info
    pathway_to_genes = get_pathway_info(nodes=None)
    print pathway_to_genes.keys()
    # Get drug similarity info
    combination_to_similarity = get_drug_similarity(drug_to_values)
    print combination_to_similarity.items()[:10]
    #values = combination_to_similarity.items()
    #values.sort(key=lambda x: x[1])
    #print values[-20:]
    task = CONFIG.get("task")
    if task.endswith("-train"):
	out_file = CONFIG.get("feature_file_train")
    elif task == "ch1-test":
	out_file = CONFIG.get("feature_file_test_ch1")
    elif task == "ch2-test":
	out_file = CONFIG.get("feature_file_test_ch2")
    else:
	raise ValueError("Uknown task: " + task)
    f = open(out_file, 'w')
    features = ["gexpA", "gexpB", "mutA", "mutB", "cnvA", "cnvB", "sim.target", "sim.chemical", "guild.med", "guild.max", "kegg.in", "kegg.gexp.med", "kegg.gexp.max", "kegg.mut.med", "kegg.mut.max", "kegg.cnv.med", "kegg.cnv.max", "kegg.cnvA", "kegg.cnvB", "cosmic.in", "cosmic.gexp.med", "cosmic.gexp.max", "cosmic.mut.med", "cosmic.mut.max", "cosmic.cnv.med", "cosmic.cnv.max", "cosmic.cnvA", "cosmic.cnvB"]
    f.write("comb.id cell.line %s\n" % " ".join(features))
    drugs = drug_to_values.keys()
    seen_combinations = set()
    for i, drug1 in enumerate(drugs):
	for j, drug2 in enumerate(drugs):
	    if i >= j:
		continue
	    #comb_id = ".".join(sorted([drug1, drug2]))
	    if drug1.lower() < drug2.lower():
		comb_id = "%s.%s" % (drug1, drug2)
	    else:
		comb_id = "%s.%s" % (drug2, drug1)
	    if comb_id not in combination_to_values:
		continue
	    if task != "ch2-test": # task.startswith("ch1-") or task.endswith("-train"):
		cell_line_to_mono_values = combination_to_values[comb_id]
	    targets1 = drug_to_values[drug1][0]
	    targets2 = drug_to_values[drug2][0]
	    for cell_line in cell_line_to_value:
		if task != "ch2-test":
		    if cell_line not in cell_line_to_mono_values:
			continue
		print comb_id, cell_line, drug1, drug2
		seen_combinations.add((comb_id, cell_line))
		feature_values = []
		# GEXP
		for targets in (targets1, targets2):
		    indices = []
		    for target in targets:
			if target in gene_to_idx:
			    indices.append(gene_to_idx[target])
		    if len(indices) == 0 or cell_line not in cell_line_to_idx:
			val = "NA"
		    else:
			val = numpy.median(gexp_norm[indices, cell_line_to_idx[cell_line]])
		    feature_values.append(val)
		# MUT
		for targets in (targets1, targets2):
		    values = []
		    for target in targets:
			if target in gene_to_cell_line_to_mutation:
			    d = gene_to_cell_line_to_mutation[target]
			    if cell_line in d: 
				values.append(d[cell_line])
		    if len(values) == 0:
			val = "NA"
		    else:
			val = numpy.max(values)
		    feature_values.append(val)
		# CNV
		for targets in (targets1, targets2):
		    values = []
		    for target in targets:
			if target in gene_to_cell_line_to_cnv:
			    d = gene_to_cell_line_to_cnv[target]
			    if cell_line in d: 
				values.append(d[cell_line])
		    if len(values) == 0:
			val = "NA"
		    else:
			val = numpy.max(values)
		    feature_values.append(val)
		# SIMILARITY
		if comb_id in combination_to_similarity:
		    vals = combination_to_similarity[comb_id]
		else:
		    vals = ["NA"] * 2
		feature_values.extend(vals)
		# GUILD
		if comb_id not in combination_to_guild_values:
		    vals = ["NA"] * 2
		else:
		    values_guild = combination_to_guild_values[comb_id]
		    vals = [ numpy.median(values_guild), numpy.max(values_guild) ]
		feature_values.extend(vals)
		# KEGG / COSMIC
		for genes in (pathway_to_genes["kegg"], pathway_to_genes["census"]):
		    # INVOLVEMENT
		    val = 0
		    for targets in (targets1, targets2):
			val += int(len(targets & genes) > 0) 
		    feature_values.append(val)
		    # GEXP
		    for target in genes:
			if target in gene_to_idx:
			    indices.append(gene_to_idx[target])
		    if len(indices) == 0 or cell_line not in cell_line_to_idx:
			vals = ["NA"] * 2
		    else:
			values = gexp_norm[indices, cell_line_to_idx[cell_line]]
			val = [ numpy.median(values), numpy.max(values) ]
		    feature_values.extend(vals)
		    # MUT
		    values = []
		    for target in genes:
			if target in gene_to_cell_line_to_mutation:
			    d = gene_to_cell_line_to_mutation[target]
			    if cell_line in d: 
				values.append(d[cell_line])
		    if len(values) == 0:
			vals = ["NA"] * 2
		    else:
			vals = [ numpy.median(values), numpy.max(values) ]
		    feature_values.extend(vals)
		    # CNV
		    values = []
		    for target in genes:
			if target in gene_to_cell_line_to_cnv:
			    d = gene_to_cell_line_to_cnv[target]
			    if cell_line in d: 
				values.append(d[cell_line])
		    if len(values) == 0:
		        vals = ["NA"] * 2
		    else:
		        vals = [ numpy.median(values), numpy.max(values) ]
		    feature_values.extend(vals)
		    # CNV target
		    for targets in (targets1, targets2):
			values = []
			for target in targets & genes:
			    if target in gene_to_cell_line_to_cnv:
				d = gene_to_cell_line_to_cnv[target]
				if cell_line in d: 
				    values.append(d[cell_line])
			if len(values) == 0:
			    val = "NA"
			else:
			    val = numpy.mean(values)
			feature_values.append(val)
		f.write("%s %s %s\n" % (comb_id, cell_line, " ".join(map(str, feature_values))))
    f.close()
    for comb_id, cell_line_to_mono_values in combination_to_values.iteritems():
	for cell_line in cell_line_to_mono_values:
	    if (comb_id, cell_line) not in seen_combinations:
		print comb_id, cell_line
    return


def get_guild_based_synergy():
    # Get network
    network = wrappers.get_network(network_file = CONFIG.get("network_file"), only_lcc = True)
    nodes = set(network.nodes())
    #create_edge_file(network)
    # Get drug info
    drug_to_values = get_drug_info(nodes=nodes)
    #print drug_to_values.keys()
    # Check individual drugs
    #guild_drugs(drug_to_values, nodes)
    # Check pairwise drug combinations
    #guild_combinations(drug_to_values, nodes)
    # Get synergy info
    combination_to_values = get_synergy_info()
    # Get gexp info
    gexp_norm, gene_to_idx, cell_line_to_idx = wrappers.get_expression_info(gexp_file = CONFIG.get("gexp_file"), process=set(["z", "abs"]), dump_file = CONFIG.get("gexp_dump"))
    # Check synergy between known pairs
    combination_to_guild_values = get_guild_based_synergy_scores(drug_to_values.keys(), None, None, None) #gexp_norm, gene_to_idx, cell_line_to_idx)
    out_file = CONFIG.get("guild_file")
    f = open(out_file, 'w')
    f.write("comb.id cell.line max.a max.b med mean sd max min syn\n")
    for comb_id, cell_line_to_values in combination_to_values.iteritems():
	#drug1, drug2 = comb_id.split(".")
	for cell_line, vals in cell_line_to_values.iteritems():
	    values = combination_to_guild_values[comb_id]
	    max_a, max_b, synergy = vals
	    f.write("%s %s %f %f %f %f %f %f %f %s\n" % (comb_id, cell_line, max_a, max_b, numpy.median(values), numpy.mean(values), numpy.std(values), numpy.max(values), numpy.min(values), synergy))
    f.close()
    return


def get_pathway_info(nodes):
    pathway_to_genes = wrappers.get_pathway_info(pathway_file = CONFIG.get("pathway_file"), prefix = CONFIG.get("pathway_source"), nodes = nodes)
    print len(pathway_to_genes) 
    #values = map(lambda x: (x[0], len(x[1])), pathway_to_genes.items())
    #values.sort(key=lambda x: x[1])
    #print values[-10:]
    #genes = pathway_to_genes["kegg_pathways_in_cancer"]
    #genes_merged = set()
    pathways_to_include = ["pathways_in_cancer", "aminoacyl-tRNA biosynthesis", "MAPK signaling pathway", "NF-kappa B signaling pathway"]
    pathways_to_include += ["Cell Cycle", "p53 signaling pathway", "Apoptosis", "TGF-beta signaling pathway"]
    #pathways_to_include += ["colorectal cancer", "small cell lung cancer", "non small cell lung cancer", "prostate cancer"]
    pathways_to_include = [ pathway.lower().replace("-", "_").replace(" ", "_") for pathway in pathways_to_include ]
    pathway_to_geneids_mod = {}
    for key, values in pathway_to_genes.iteritems():
	if key.find("cancer") != -1 and key != "kegg_pathways_in_cancer":
	    print key, len(values)
	    #genes_merged |= pathway_to_genes[key]
	elif key.find("growth") != -1 or key.find("apoptosis") != -1:
	    print key, len(values)
    genes = set()
    for pathway in pathways_to_include:
	try:
	    #pathway_to_geneids_mod[pathway] = pathway_to_genes["kegg_" + pathway]
	    genes |= pathway_to_genes["kegg_" + pathway]
	except:
	    print "Pathway not found!", pathway
    pathway_to_geneids_mod["kegg"] = genes
    conditions_to_include = ["census"] #"lung", "breast", "colo", "prostate", "gastrointestinal"])
    condition_to_genes = get_census_gene_info()
    for pathway in conditions_to_include:
	pathway_to_geneids_mod[pathway] = condition_to_genes[pathway]
    #print len(genes), len(genes_merged), len(genes & genes_merged)
    return pathway_to_geneids_mod


def get_census_gene_info():
    file_name = CONFIG.get("census_file")
    condition_to_genes = {}
    f = open(file_name)
    reader = csv.reader(f, delimiter=",", quotechar='"')
    header = reader.next()
    #print header
    header_to_idx = dict((word, i) for i, word in enumerate(header))
    conditions_to_include = set(["lung", "breast", "colo", "prostate", "gastrointestinal"])
    #print header_to_idx
    genes = set()
    for row in reader:
	gene = row[header_to_idx['Gene Symbol']]
	type_somatic = row[header_to_idx['Tumour Types(Somatic)']]
	type_germline = row[header_to_idx['Tumour Types(Germline)']]
	#if any(map(conditions_to_include, lambda x: type_somatic.find(x) != -1)) or any(map(conditions_to_include, lambda x: type_somatic.find(x) != -1)):
	for condition in conditions_to_include:
	    if type_somatic.find(condition) != -1 or type_germline.find(condition) != -1:
		condition_to_genes.setdefault(condition, set()).add(gene)
	genes.add(gene)
    condition_to_genes["census"] = genes
    f.close()
    return condition_to_genes


def get_cell_line_info():
    file_name = CONFIG.get("cell_line_file")
    cell_line_to_value = {}
    f = open(file_name)
    reader = csv.reader(f, delimiter=",", quotechar='"')
    header = reader.next()
    header_to_idx = dict((word, i) for i, word in enumerate(header))
    for row in reader:
	cell = row[header_to_idx['Sanger.Name']]
	tissue = row[header_to_idx['Tissue..General.']]
	cell_line_to_value[cell] = tissue
    f.close()
    return cell_line_to_value
	

def get_mutation_info():
    dump_file = CONFIG.get("mutation_dump")
    if os.path.exists(dump_file):
	gene_to_cell_line_to_value = cPickle.load(open(dump_file))
	return gene_to_cell_line_to_value 
    file_name = CONFIG.get("mutation_file")
    gene_to_cell_line_to_value = {}
    f = open(file_name)
    reader = csv.reader(f, delimiter=",", quotechar='"')
    header = reader.next()
    #print header
    header_to_idx = dict((word, i) for i, word in enumerate(header))
    #print header_to_idx
    for row in reader:
	gene = row[header_to_idx['Gene.name']]
	gene = gene.split(".")[0]
	gene = gene.split("_")[0]
	cell_line_to_value = gene_to_cell_line_to_value.setdefault(gene, {})
	cell_line = row[header_to_idx['cell_line_name']]
	description = row[header_to_idx['Mutation.Description']]
	prediction = row[header_to_idx['FATHMM.prediction']]
	mut = None
	if description == "Substitution - coding silent" or description == "Unknown":
	    mut = 0
	elif prediction =="CANCER":
	    mut = 2
	else:
	    mut = 1
	cell_line_to_value[cell_line] = mut
    f.close()
    cPickle.dump(gene_to_cell_line_to_value, open(dump_file, 'w')) 
    return gene_to_cell_line_to_value


def get_cnv_info(file_name):
    dump_file = CONFIG.get("cnv_dump")
    if os.path.exists(dump_file):
	gene_to_cell_line_to_value = cPickle.load(open(dump_file))
	return gene_to_cell_line_to_value 
    gene_to_cell_line_to_value = {}
    f = open(file_name)
    reader = csv.reader(f, delimiter=",", quotechar='"')
    header = reader.next()
    #print header
    header_to_idx = dict((word, i) for i, word in enumerate(header))
    #print header_to_idx
    for row in reader:
	ch = row[header_to_idx['chr_GRCh37']]
	if ch == "Y":
	    continue
	gene = row[header_to_idx['gene']]
	gene = gene.split(".")[0]
	gene = gene.split("_")[0]
	cell_line_to_value = gene_to_cell_line_to_value.setdefault(gene, {})
	cell_line = row[header_to_idx['cell_line_name']]
	cnv_min = int(row[header_to_idx['min_cn_GRCh37']])
	cnv_max = int(row[header_to_idx['max_cn_GRCh37']])
	cnv = (cnv_min + cnv_max) / 2.0
	cell_line_to_value[cell_line] = cnv
    f.close()
    cPickle.dump(gene_to_cell_line_to_value, open(dump_file, 'w')) 
    return gene_to_cell_line_to_value


def get_drug_info(nodes=None):
    drug_file = CONFIG.get("drug_file")
    drug_to_values = {}
    f = open(drug_file)
    reader = csv.reader(f, delimiter=',', quotechar='"')
    header = reader.next()
    for row in reader:
	drug = row[0]
	targets = set(row[1].split(","))
	if nodes is not None:
	    targets &= nodes
	#if len(targets) == 0:
	#    continue
	smiles = row[-2]
	smiles = list(sorted(map(lambda x: (len(x), x), smiles.split(";"))))[-1][1]
	drug_to_values[drug] = (targets, smiles)
	#print drug, targets
    f.close()
    return drug_to_values


def get_drug_similarity(drug_to_values):
    combination_to_similarity = {}
    drugs = drug_to_values.keys()
    for i, drug1 in enumerate(drugs):
	for j, drug2 in enumerate(drugs):
	    if i >= j:
		continue
	    comb_id = ".".join(sorted([drug1, drug2]))
	    targets1, smiles1 = drug_to_values[drug1]
	    targets2, smiles2 = drug_to_values[drug2]
	    targets_common = targets1 & targets2
	    values = []
	    for method in ("target", "chemical"):
		if method == "target":
		    if len(targets1) == 0 or len(targets2) == 0:
			#continue
			d = "NA" 
		    else:
			d = len(targets_common) / float(len(targets1|targets2)) # float(min(len(targets1), len(targets2)))
		    #values.append(d)
		    combination_to_similarity.setdefault(comb_id, []).append(d)
		elif method == "chemical":
		    if len(smiles1) == 0 or len(smiles2) == 0:
			#continue
			d = "NA"
		    else:
			ind = indigo.Indigo()
			m = ind.loadMolecule(smiles1) 
			m.aromatize()
			fp = m.fingerprint("sim") # sub
			m2 = ind.loadMolecule(smiles2) 
			m2.aromatize() # Aromatize molecules in case they are not in aromatic form
			fp2 = m2.fingerprint("sim") # Calculate similarity between "similarity" fingerprints
			d = ind.similarity(fp, fp2, "tanimoto") 
		    #values.append(d/2) # Weight chemical similarity less (half of target similiarity)
		    combination_to_similarity.setdefault(comb_id, []).append(d)
		    #print group, group2, "Tanimoto: %s" % (d) 
		else:
		    raise ValueError("Uknown method: %s" % method)
	    #combination_to_similarity[comb_id] = numpy.mean(values)
    return combination_to_similarity


def get_synergy_info():
    task = CONFIG.get("task")
    if task.endswith("-train"):
	combination_file = CONFIG.get("combination_file_train")
    elif task == "ch1-test":
	combination_file = CONFIG.get("combination_file_test_ch1")
    elif task == "ch2-test":
	combination_file = CONFIG.get("combination_file_test_ch2")
    else:
	raise ValueError("Uknown task: " + task)
    combination_to_values = {}
    f = open(combination_file )
    reader = csv.reader(f, delimiter=',', quotechar='"')
    header = reader.next()
    for row in reader:
	cell_line, drug1, drug2 = row[:3]
	max_a, max_b = map(float, row[3:5])
	#"IC50_A","H_A","Einf_A","IC50_B","H_B","Einf_B"
	synergy, qa, comb_id = row[-3:]
	qa = int(qa)
	if qa != 1:
	    continue
	#synergy = float(synergy)
	if comb_id != "%s.%s" % (drug1, drug2):
	    print comb_id, drug1, drug2
	cell_line_to_values = combination_to_values.setdefault(comb_id, {})
	if cell_line in cell_line_to_values: 
	    print "Overwriting", comb_id, cell_line, max_a, max_b
	cell_line_to_values[cell_line] = (max_a, max_b, synergy)
    f.close()
    return combination_to_values 


def create_edge_file(network):
    network_lcc_file = CONFIG.get("network_file") + ".lcc"
    f = open(network_lcc_file, 'w')
    for u,v in network.edges():
	f.write("%s 1 %s\n" % (u, v))
    f.close()
    return


def create_node_file(drug, targets, nodes):
    output_dir = CONFIG.get("output_dir") + "/"
    node_file = "%s%s.node" % (output_dir, drug)
    f = open(node_file, 'w')
    for node in nodes:
	if node in targets:
	    score = 1
	else:
	    score = 0.01
	f.write("%s %f\n" % (node, score))
    f.close()
    return


def get_guild_based_synergy_scores(drugs, gexp_norm, gene_to_idx, cell_line_to_idx):
    dump_file = CONFIG.get("guild_dump")
    if os.path.exists(dump_file):
	combination_to_values = cPickle.load(open(dump_file))
	return combination_to_values
    output_dir = CONFIG.get("output_dir") + "/"
    use_expression = int(CONFIG.get("use_expression"))
    combination_to_values = {}
    for i, drug1 in enumerate(drugs):
	for j, drug2 in enumerate(drugs):
	    if i >= j:
		continue
	    comb_id = ".".join(sorted([drug1, drug2]))
	    #print comb_id, drug1, drug2
	    drug1_file = "%s%s.ns" % (output_dir, drug1)
	    drug2_file = "%s%s.ns" % (output_dir, drug2)
	    comb_file = "%s%s.%s.ns" % (output_dir, drug1, drug2)
	    if not os.path.exists(comb_file):
		comb_file = "%s%s.%s.ns" % (output_dir, drug2, drug1)
		if not os.path.exists(comb_file):
		    #print "Combination not found!"
		    continue
	    drug1_to_score, values1 = get_top_genes_and_scores(drug1_file)
	    drug2_to_score, values2 = get_top_genes_and_scores(drug2_file)
	    genes_common = set(zip(*values1)[0]) & set(zip(*values2)[0])
	    #print comb_id, drug1, drug2, len(genes_common)
	    if len(genes_common) == 0:
		continue
	    comb_to_score = guild_utilities.get_node_to_score(comb_file)
	    values = []
	    if use_expression == 0:
		for gene in genes_common:
		    val = comb_to_score[gene] - (drug1_to_score[gene] + drug2_to_score[gene]) / 2
		    values.append(val)
	    elif use_expression == 1:
		genes_common_new = set()
		for gene in genes_common:
		    if gene not in gene_to_idx or cell_line not in cell_line_to_idx:
			continue
		    exp = gexp_norm[gene_to_idx[gene], cell_line_to_idx[cell_line]]
		    if exp < 1:
			continue
		    genes_common_new.add(gene)
		genes_common = genes_common_new
		for gene in genes_common:
		    val = comb_to_score[gene] - (drug1_to_score[gene] + drug2_to_score[gene]) / 2
		    values.append(val)
	    elif use_expression == 2:
		for gene in genes_common:
		    if gene not in gene_to_idx or cell_line not in cell_line_to_idx:
			continue
		    exp = gexp_norm[gene_to_idx[gene], cell_line_to_idx[cell_line]]
		    values.append(exp)
	    if len(values) == 0:
		continue
	    combination_to_values[comb_id] = values
    cPickle.dump(combination_to_values, open(dump_file, 'w')) 
    return combination_to_values 


def get_top_genes_and_scores(score_file):
    node_to_score = guild_utilities.get_node_to_score(score_file)
    #print node_to_score.items()[:3]
    values = node_to_score.items()
    values.sort(key=lambda x: x[1])
    return node_to_score, values[-500:]


def guild_combinations(drug_to_values, nodes):
    drugs = drug_to_values.keys()
    flag = False # True #
    for i, drug1 in enumerate(drugs):
	if drug1 == "NAE":
	    flag = True
	if drug1 == "Docetaxel":
	    flag = False
	if not flag:
	    continue
	for j, drug2 in enumerate(drugs):
	    if i < j:
		targets1, smiles1 = drug_to_values[drug1]
		targets2, smiles2 = drug_to_values[drug2]
		targets = targets1 | targets2
		if len(targets) == max(len(targets1), len(targets2)):
		    #print "Combination redundant!"
		    continue
		combination = "%s.%s" % (drug1, drug2)
		#print combination, targets
		#if len(targets & nodes) == 0:
		#   print "Not in network!"
		#   continue
		run_guild(combination, targets, qname="all.q")
    return


def guild_drugs(drug_to_values, nodes):
    for drug, values in drug_to_values.iteritems():
	targets = values[0]
	print drug, targets
	#if len(targets & nodes) == 0:
	#   print "Not in network!"
	#   continue
	run_guild(drug, targets, qname="all.q")
    return


def run_guild(drug, targets, qname=None):
    # Create node file
    create_node_file(drug, targets, nodes)
    # Get parameters
    executable_path = CONFIG.get("guild_path") 
    output_dir = CONFIG.get("output_dir") + "/"
    network_lcc_file = CONFIG.get("network_file") + ".lcc"
    output_file = "%s%s.ns" % (output_dir, drug)
    node_file = "%s%s.node" % (output_dir, drug)
    n_repetition = int(CONFIG.get("n_repetition"))
    n_iteration = int(CONFIG.get("n_iteration"))
    # Get and run the GUILD command
    score_command = executable_path + ' -s s -n "%s" -e "%s" -o "%s" -r %d -i %d' % (node_file, network_lcc_file, output_file, n_repetition, n_iteration)
    #print strftime("%H:%M:%S - %d %b %Y") #, score_command
    if qname is None:
	os.system(score_command)
    else:
	#os.system("qsub -cwd -o out -e err -q %s -N %s -b y %s" % (qname, scoring_type, score_command))
	print "qsub -cwd -o out -e err -q %s -N guild_%s -b y %s" % (qname, drug, score_command)
    return


if __name__ == "__main__":
    main()
    
