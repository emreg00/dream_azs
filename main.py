
from toolbox import configuration, guild_utilities, wrappers, network_utilities
from rdbox import drug_info
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
    # Get SEA & STITCH predictions
    #print get_drug_info_with_predicted_targets()
    return


def get_drug_info_with_predicted_targets():
    file_name = CONFIG.get("drug_file") + ".predicted"
    # Get drug info
    drug_to_values = get_drug_info(nodes=None)
    if os.path.exists(file_name):
	f = open(file_name)
	f.readline()
	drug_to_targets_predicted = {}
	for line in f:
	    drug, targets, targets_org = line.strip("\n").split("\t")
	    if targets != "":
		drug_to_targets_predicted[drug] = set(targets.split(", "))
    else:
	#smiles_list = ["[H][C@]12C[C@@H](O)C=C[C@]11CCN(C)CC3=C1C(O2)=C(OC)C=C3"]
	smiles_list = []
	smiles_to_drug = {}
	for drug, values in drug_to_values.iteritems():
	    smiles = values[1] 
	    if smiles.strip() == "":
		continue
	    smiles_list.append(smiles)
	    smiles_to_drug[smiles] = drug
	# Get SMILES similarity based predictions
	#smiles_to_geneids = drug_info.get_chemical_similarity_based_target_predictions(CONFIG, smiles_list, cutoff=0.7, method="any_smiles") # majority_above
	# Fisher's based test for taking into account of abundance of same target in many similar / dissimilar smiles
	smiles_to_geneids = drug_info.get_chemical_similarity_based_target_predictions(CONFIG, smiles_list, cutoff=0.6, method="fishers") # 0.7 
	#print len(smiles_to_geneids), smiles_to_geneids.items()[:5]
	# Get gene names
	geneid_to_names, name_to_geneid = wrappers.get_geneid_symbol_mapping(CONFIG.get("id_mapping_file"))
	f = open(file_name, 'w')
	f.write("Drug\tTargets (Predicted)\tTargets (Original)\n")
	for smiles in smiles_list:
	    genes = set()
	    if smiles in smiles_to_geneids:
		for geneid in smiles_to_geneids[smiles]:
		    genes |= geneid_to_names[geneid]
	    genes = list(genes)
	    genes.sort()
	    #f.write("%s\t%s\t%s\n" % (smiles_to_drug[smiles], smiles, ", ".join(genes)))
	    targets_org = []
	    if smiles in smiles_to_drug:
		drug = smiles_to_drug[smiles]
		if drug in drug_to_values:
		    targets_org = drug_to_values[drug][0]
	    f.write("%s\t%s\t%s\n" % (smiles_to_drug[smiles], ", ".join(genes), ", ".join(targets_org)))
	f.close()
    # Can check targets from STITCH for benchmarking
    for drug, values in drug_to_values.iteritems():
	if drug in drug_to_targets_predicted:
	    targets_predicted = drug_to_targets_predicted[drug]
	    values[0] = values[0] | targets_predicted
    return drug_to_values


def create_feature_file():
    # Get drug info
    drug_to_values = get_drug_info(nodes=None) 
    #drug_to_values = get_drug_info_with_predicted_targets() # worsens the performance
    #for drug, values in drug_to_values.iteritems():
    #	if values[1] != "":
    #	    print drug, values[1]
    #print drug_to_values.items()[:3]
    # Get cell line info
    cell_line_to_value = get_cell_line_info()
    #print cell_line_to_value.items()[:3]
    # Get synergy info
    combination_to_values = get_synergy_info()
    cell_line_to_synergy, combination_to_synergy = get_synergy_values_per_cell_line_and_combination()
    #print combination_to_values.items()[:3]
    # Get gexp info
    gexp_norm, gene_to_idx, cell_line_to_idx = wrappers.get_expression_info(gexp_file = CONFIG.get("gexp_file"), process=set(["z"]), dump_file = CONFIG.get("gexp_dump")) # process=set(["z", "abs"])
    #values = [gene_to_idx["TSPAN6"], gene_to_idx["TNMD"]] 
    #print gexp_norm[values, cell_line_to_idx["647-V"]]
    #print "TNMD @ 647-V", gexp_norm[gene_to_idx["TNMD"], cell_line_to_idx["647-V"]] 
    # Get methylation info
    meth, meth_gene_to_idx, meth_cell_line_to_idx = wrappers.get_expression_info(gexp_file = CONFIG.get("methylation_file"), process=set(["z"]), dump_file = CONFIG.get("methylation_dump")) 
    #print "A1BG @ 647-V", meth[meth_gene_to_idx["A1BG"], meth_cell_line_to_idx["647-V"]]
    # Get mutation info
    gene_to_cell_line_to_mutation = get_mutation_info() 
    #print gene_to_cell_line_to_mutation.items()[:3]
    # Get CNV info
    #gene_to_cell_line_to_cnv = {}
    gene_to_cell_line_to_cnv = get_cnv_info(CONFIG.get("cnv_file")) 
    #print gene_to_cell_line_to_cnv.items()[:2]
    # Get cancer gene & pathway info
    pathway_to_genes = get_pathway_info(nodes=None)
    #print pathway_to_genes.keys()
    genes_pathway = set()
    for genes in (pathway_to_genes["kegg"], pathway_to_genes["census"]):
	genes_pathway |= genes
    genes_pathway = list(genes_pathway)
    # Get network
    network = wrappers.get_network(network_file = CONFIG.get("network_file"), only_lcc = True)
    network_nodes = set(network.nodes())
    # Get GUILD info
    #combination_to_guild_values = {}
    combination_to_guild_values = get_guild_based_synergy_scores(drug_to_values.keys(), genes_pathway, gexp_norm, gene_to_idx, cell_line_to_idx)
    #print combination_to_guild_values.items()[:2]
    # Get drug similarity info
    combination_to_similarity = get_drug_similarity(drug_to_values)
    #print combination_to_similarity.items()[:5]
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
    features = ["gexpA.med", "gexpA.amed", "gexpB.med", "gexpB.amed", "mutA", "mutB", "cnvA", "cnvB", "metA.med", "metA.amed", "metB.med", "metB.amed", "cell.med", "comb.med", "sim.target", "sim.chemical", "kA", "kB", "dAB", "guild.common", "guild.med", "guild.max", "kegg.inA", "kegg.inB", "kegg.gexp.med", "kegg.gexp.max", "kegg.mut.med", "kegg.mut.max", "kegg.cnv.med", "kegg.cnv.max", "kegg.cnvA", "kegg.cnvB", "cosmic.inA", "cosmic.inB", "cosmic.gexp.med", "cosmic.gexp.max", "cosmic.mut.med", "cosmic.mut.max", "cosmic.cnv.med", "cosmic.cnv.max", "cosmic.cnvA", "cosmic.cnvB"]
    drugs = drug_to_values.keys()
    seen_combinations = set()
    # Get all targets
    targets_all = set()
    for i, drug1 in enumerate(drugs):
	for j, drug2 in enumerate(drugs):
	    if i >= j:
		continue
	    targets1 = drug_to_values[drug1][0]
	    targets2 = drug_to_values[drug2][0]
	    targets_all |= targets1 | targets2
    targets_all = list(targets_all)
    target_to_idx = dict((target, i) for i, target in enumerate(targets_all))
    # Get all pathway genes
    #targets_pathway = set()
    #for genes in (pathway_to_genes["kegg"], pathway_to_genes["census"]): 
    #	targets_pathway |= genes
    #targets_pathway = list(targets_pathway & set(targets_all))
    #pathway_target_to_idx = dict((gene, i) for i, gene in enumerate(targets_pathway))
    targets_pathway = []
    for pathway, genes in pathway_to_genes.iteritems():
	targets_pathway.append(pathway)
    # Create header
    features = map(lambda x: ".g" + x, targets_all) + map(lambda x: ".m" + x, targets_all) + map(lambda x: ".c" + x, targets_all) + map(lambda x: ".z" + x, targets_all) + map(lambda x: ".e" + x, targets_all) + map(lambda x: ".p" + x.replace("_", "."), targets_pathway) + features
    f.write("comb.id cell.line %s\n" % " ".join(features))
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
		#print comb_id, cell_line, drug1, drug2
		seen_combinations.add((comb_id, cell_line))
		feature_values = []
		# GEXP categorized
		values = [ 0 ] * len(targets_all)
		if cell_line in cell_line_to_idx:
		    for targets in (targets1, targets2):
			for target in targets:
			    if target in gene_to_idx:
				values[target_to_idx[target]] += gexp_norm[gene_to_idx[target], cell_line_to_idx[cell_line]]
		feature_values.extend(values)
		#print len(feature_values) 
		# MUT categorized
		values = [ 0 ] * len(targets_all)
		for targets in (targets1, targets2):
		    for target in targets:
			if target in gene_to_cell_line_to_mutation:
			    d = gene_to_cell_line_to_mutation[target]
			    if cell_line in d: 
				values[target_to_idx[target]] += d[cell_line]
		feature_values.extend(values)
		#print len(feature_values)
		# CNV & ZYG categorized
		values = [ 0 ] * len(targets_all)
		values2 = [ 0 ] * len(targets_all)
		for targets in (targets1, targets2):
		    for target in targets:
			if target in gene_to_cell_line_to_cnv:
			    d = gene_to_cell_line_to_cnv[target] 
			    if cell_line in d: 
				# d contains cnv value and zygosity
				values[target_to_idx[target]] += d[cell_line][0]
				values2[target_to_idx[target]] += d[cell_line][1]
		feature_values.extend(values) 
		feature_values.extend(values2) 
		# METH categorized
		values = [ 0 ] * len(targets_all)
		if cell_line in meth_cell_line_to_idx:
		    for targets in (targets1, targets2):
			for target in targets:
			    if target in meth_gene_to_idx:
				values[target_to_idx[target]] += meth[meth_gene_to_idx[target], meth_cell_line_to_idx[cell_line]]
		feature_values.extend(values)	
		#print len(feature_values)
		# KEGG / COSMIC INVOLVEMENT categorized
		values = [ 0 ] * len(targets_pathway)
		#for targets in (targets1, targets2):
		    #for target in targets:
			#if target in pathway_target_to_idx:
			#    values[pathway_target_to_idx[target]] += 1
		for k, pathway in enumerate(targets_pathway):
		    values[k] += len(targets1 & pathway_to_genes[pathway])
		    values[k] += len(targets2 & pathway_to_genes[pathway])
		feature_values.extend(values) 
		#print len(feature_values) 
		# GEXP
		for targets in (targets1, targets2):
		    indices = []
		    for target in targets:
			if target in gene_to_idx:
			    indices.append(gene_to_idx[target])
		    if len(indices) == 0 or cell_line not in cell_line_to_idx:
			vals = ["NA"] * 2
		    else:
			values = gexp_norm[indices, cell_line_to_idx[cell_line]]
			vals = [numpy.median(values), numpy.median(numpy.abs(values))] #values.flat[numpy.abs(values).argmax()]]
		    feature_values.extend(vals)
		    #print len(feature_values)
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
		    #print len(feature_values) 
		# CNV
		for targets in (targets1, targets2):
		    values = []
		    for target in targets:
			if target in gene_to_cell_line_to_cnv:
			    d = gene_to_cell_line_to_cnv[target]
			    if cell_line in d: 
				values.append(d[cell_line][0])
		    if len(values) == 0:
			val = "NA"
		    else:
			val = numpy.max(values)
		    feature_values.append(val)
		    #print len(feature_values)
		# METH
		for targets in (targets1, targets2):
		    indices = []
		    for target in targets:
			if target in meth_gene_to_idx:
			    indices.append(meth_gene_to_idx[target])
		    if len(indices) == 0 or cell_line not in meth_cell_line_to_idx:
			vals = ["NA"] * 2
		    else:
			#val = numpy.median(numpy.abs(meth[indices, meth_cell_line_to_idx[cell_line]]))
			values = meth[indices, meth_cell_line_to_idx[cell_line]]
			vals = [numpy.median(values), numpy.median(numpy.abs(values))] 
		    feature_values.extend(vals)
		    #print len(feature_values)
		# MEDIAN SYNERGY per cell line / combination
		vals = ["NA"] * 2
		if cell_line in cell_line_to_synergy:
		    vals[0] = cell_line_to_synergy[cell_line]
		if comb_id in combination_to_synergy:
		    vals[1] = combination_to_synergy[comb_id]
		# SIMILARITY
		vals = ["NA"] * 2
		if comb_id in combination_to_similarity:
		    vals = combination_to_similarity[comb_id]
		feature_values.extend(vals)
		#print len(feature_values) 
		# Interaction network based features (degree A/B)
		for targets in (targets1, targets2):
		    values = []
		    for target in targets:
			if target in network_nodes:
			    d = network.degree(target)
			    values.append(d)
		    if len(values) == 0:
			val = "NA"
		    else:
			val = numpy.max(values)
		    feature_values.append(val)
		# Interaction network based distance between A-B
		values = []
		for target1 in targets1:
		    if target1 not in network_nodes:
			continue
		    for target2 in targets2:
			if target2 not in network_nodes:
			    continue
			d = network_utilities.get_shortest_path_length_between(network, target1, target2)
			values.append(d)
		if len(values) == 0:
		    val = "NA"
		else:
		    val = numpy.min(values)
		feature_values.append(val)
		# GUILD
		vals = ["NA"] * 3
		if comb_id in combination_to_guild_values:
		    values_guild = combination_to_guild_values[comb_id]
		    if cell_line in values_guild:
			vals = values_guild[cell_line]
		feature_values.extend(vals)
		#print len(feature_values)
		# KEGG / COSMIC
		for genes in (pathway_to_genes["kegg"], pathway_to_genes["census"]):
		    # INVOLVEMENT
		    val = 0
		    for targets in (targets1, targets2):
			val = len(targets & genes)
			feature_values.append(val)
		    #print len(feature_values) 
		    # GEXP
		    for target in genes:
			if target in gene_to_idx:
			    indices.append(gene_to_idx[target])
		    if len(indices) == 0 or cell_line not in cell_line_to_idx:
			vals = ["NA"] * 2
		    else:
			values = gexp_norm[indices, cell_line_to_idx[cell_line]]
			vals = [ numpy.median(values), numpy.max(values) ]
		    feature_values.extend(vals)
		    #print len(feature_values) 
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
		    #print len(feature_values) 
		    # CNV
		    values = []
		    for target in genes:
			if target in gene_to_cell_line_to_cnv:
			    d = gene_to_cell_line_to_cnv[target]
			    if cell_line in d: 
				values.append(d[cell_line][0])
		    if len(values) == 0:
		        vals = ["NA"] * 2
		    else:
		        vals = [ numpy.median(values), numpy.max(values) ]
		    feature_values.extend(vals)
		    #print len(feature_values) 
		    # CNV target
		    for targets in (targets1, targets2):
			values = []
			for target in targets & genes:
			    if target in gene_to_cell_line_to_cnv:
				d = gene_to_cell_line_to_cnv[target]
				if cell_line in d: 
				    values.append(d[cell_line][0])
			if len(values) == 0:
			    val = "NA"
			else:
			    val = numpy.mean(values)
			feature_values.append(val)
			#print len(feature_values) 
		f.write("%s %s %s\n" % (comb_id, cell_line, " ".join(map(str, feature_values))))
		#if comb_id == "BCL2L1.Vinorelbine" and cell_line == "NCI-H1437": 
		#    return
    f.close()
    print "Not in seen combinations:"
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
    # Get gexp info
    gexp_norm, gene_to_idx, cell_line_to_idx = None, None, None
    #gexp_norm, gene_to_idx, cell_line_to_idx = wrappers.get_expression_info(gexp_file = CONFIG.get("gexp_file"), process=set(["z", "abs"]), dump_file = CONFIG.get("gexp_dump"))
    # Check individual drugs
    guild_drugs(drug_to_values, nodes, gexp_norm, gene_to_idx, cell_line_to_idx)
    # Check pairwise drug combinations
    guild_combinations(drug_to_values, nodes, gexp_norm, gene_to_idx, cell_line_to_idx)
    return 
    # Now using create_feature_file instead of below
    # Get synergy info
    combination_to_values = get_synergy_info()
    # Check synergy between known pairs
    combination_to_guild_values = get_guild_based_synergy_scores(drug_to_values.keys(), None, gexp_norm, gene_to_idx, cell_line_to_idx)
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
    #print len(pathway_to_genes) 
    #values = map(lambda x: (x[0], len(x[1])), pathway_to_genes.items())
    #values.sort(key=lambda x: x[1])
    #print values[-10:]
    #genes = pathway_to_genes["kegg_pathways_in_cancer"]
    #genes_merged = set()
    #pathways_to_include = ["pathways in cancer", "aminoacyl-tRNA biosynthesis", "MAPK signaling pathway", "NF-kappa B signaling pathway"]
    #pathways_to_include += ["Cell Cycle", "p53 signaling pathway", "Apoptosis", "TGF-beta signaling pathway"]
    pathways_to_include = CONFIG.get("pathways_to_include").split("|")
    #pathways_to_include += ["colorectal cancer", "small cell lung cancer", "non small cell lung cancer", "prostate cancer"]
    pathways_to_include = [ pathway.lower().replace("-", "_").replace(" ", "_") for pathway in pathways_to_include ]
    pathway_to_geneids_mod = {}
    for key, values in pathway_to_genes.iteritems():
	if key.find("cancer") != -1 and key != "kegg_pathways_in_cancer":
	    pass
	    #print key, len(values)
	    #genes_merged |= pathway_to_genes[key]
	elif key.find("growth") != -1 or key.find("apoptosis") != -1:
	    pass
	    #print key, len(values)
    genes = set()
    for pathway in pathways_to_include:
	try:
	    pathway_to_geneids_mod[pathway] = pathway_to_genes["kegg_" + pathway]
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
	# take into account zygosity (H:0, LOH:1, O:2)
	zygosity = row[header_to_idx['zygosity_GRCh37']]
	if zygosity == "H":
	    zygosity = 0
	elif zygosity == "L":
	    zygosity = 1
	elif zygosity == "0":
	    zygosity = 2
	else:
	    zygosity = 0
	cell_line_to_value[cell_line] = (cnv, zygosity)
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
	drug_to_values[drug] = [targets, smiles]
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
	    #comb_id = ".".join(sorted([drug1, drug2]))
	    if drug1.lower() < drug2.lower():
		comb_id = "%s.%s" % (drug1, drug2)
	    else:
		comb_id = "%s.%s" % (drug2, drug1)
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


def get_synergy_values_per_cell_line_and_combination():
    # Assigning median synergy for per cell line / drug combination
    cell_line_to_synergy = {}
    combination_to_synergy = {}
    combination_to_values = get_synergy_info(task="-train")
    for comb_id, cell_line_to_mono_values in combination_to_values.iteritems():
	for cell_line, values in cell_line_to_mono_values.iteritems():
	    cell_line_to_synergy.setdefault(cell_line, []).append(values[-1])
	    combination_to_synergy.setdefault(combination, []).append(values[-1])
    for cell_line, values in cell_line_to_synergy.items():
	cell_line_to_synergy[cell_line] = numpy.median(values)
    for combination, values in combination_to_synergy.items():
	 combination_to_synergy[combination] = numpy.median(values)
    return cell_line_to_synergy, combination_to_synergy


def get_synergy_info(task=None):
    if task is None:
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
	if task.endswith("-train"):
	    synergy = float(synergy)
	else:
	    synergy = None
	qa = int(qa)
	if qa != 1:
	    continue
	#synergy = float(synergy)
	if comb_id != "%s.%s" % (drug1, drug2):
	    print comb_id, drug1, drug2
	if drug1.lower() > drug2.lower():
	    raise ValueError("Inconsistent format!")
	cell_line_to_values = combination_to_values.setdefault(comb_id, {})
	flag = True
	if cell_line in cell_line_to_values: 
	    #print "Keeping 1/1 concentration (if exists) or max abs synergy", comb_id, cell_line, max_a, max_b, synergy, "(", cell_line_to_values[cell_line],")"
	    max_a2, max_b2, synergy2 = cell_line_to_values[cell_line]
	    if max_a2 == 1 and max_b2 == 1:
		if synergy is not None and max_a == 1 and max_b == 1 and abs(synergy) > abs(synergy2):
		    flag = True
		else:
		    flag = False
	    elif max_a == 1 and max_b == 1:
		flag = True
	    elif synergy is not None and abs(synergy) > abs(synergy2):
		flag = True
	if flag:
	    cell_line_to_values[cell_line] = [max_a, max_b, synergy]
    f.close()
    return combination_to_values 


def create_edge_file(network):
    network_lcc_file = CONFIG.get("network_file") + ".lcc"
    f = open(network_lcc_file, 'w')
    for u,v in network.edges():
	f.write("%s 1 %s\n" % (u, v))
    f.close()
    return


def create_node_file(drug, target_to_score, nodes):
    output_dir = CONFIG.get("output_dir") + "/"
    node_file = "%s%s.node" % (output_dir, drug)
    f = open(node_file, 'w')
    for node in nodes:
	if node in target_to_score:
	    score = target_to_score[node]
	else:
	    score = 0.01
	f.write("%s %f\n" % (node, score))
    f.close()
    return


def get_guild_based_synergy_scores(drugs, genes_pathway, gexp_norm, gene_to_idx, cell_line_to_idx):
    output_dir = CONFIG.get("output_dir") + "/"
    use_expression = CONFIG.get("use_expression")
    dump_file = CONFIG.get("guild_dump")
    if os.path.exists(dump_file):
	drug_to_top_genes = cPickle.load(open(dump_file))
	#return combination_to_values
    else:
	drug_to_top_genes = {}
	if use_expression != "pathway-genes":
	    for i, drug in enumerate(drugs):
		drug_file = "%s%s.ns" % (output_dir, drug)
		if not os.path.exists(drug_file): # no targets
		    drug_file = "%s%s.ns" % (output_dir, "background")
		gene_to_score, values = get_top_genes_and_scores(drug_file)
		# Top 500 genes and their scores
		for k, cell_line in enumerate(cell_line_to_idx):
		    #print i, k, drug, cell_line
		    drug_to_top_genes.setdefault(drug, {})[cell_line] = dict(values)
	else:
	    for i, drug in enumerate(drugs):
		for k, cell_line in enumerate(cell_line_to_idx):
		    #print i, k, drug, cell_line
		    drug_file = "%s%s.%s.ns" % (output_dir, drug, cell_line) 
		    gene_to_score, values = get_top_genes_and_scores(drug_file)
		    #drug_to_top_genes.setdefault(drug, {})[cell_line] = set(zip(*values)[0])
		    # Top 500 genes and their scores
		    drug_to_top_genes.setdefault(drug, {})[cell_line] = dict(values)
	cPickle.dump(drug_to_top_genes, open(dump_file, 'w')) 
    #pathway_gene_to_idx = dict((gene, i) for i, gene in enumerate(genes_pathway))
    combination_to_values = {}
    for i, drug1 in enumerate(drugs):
	#print i, drug1
	for j, drug2 in enumerate(drugs):
	    if i >= j:
		continue
	    #print j, drug2
	    #comb_id = ".".join(sorted([drug1, drug2]))
	    if drug1.lower() < drug2.lower():
		comb_id = "%s.%s" % (drug1, drug2)
	    else:
		comb_id = "%s.%s" % (drug2, drug1)
	    if use_expression != "pathway-genes":
		comb_to_score = None
		comb_file = "%s%s.%s.ns" % (output_dir, drug1, drug2)
		if os.path.exists(comb_file):
		    comb_to_score = guild_utilities.get_node_to_score(comb_file)
		else:
		    comb_file = "%s%s.%s.ns" % (output_dir, drug2, drug1)
		    if os.path.exists(comb_file):
			comb_to_score = guild_utilities.get_node_to_score(comb_file)
		    #else:
			#print "Combination not found:", comb_id
			#continue
	    cell_line_to_values = {} 
	    for cell_line in cell_line_to_idx:
		genes_to_score1 = drug_to_top_genes[drug1][cell_line]
		genes_to_score2 = drug_to_top_genes[drug2][cell_line]
		genes_common = set(genes_to_score1.keys()) & set(genes_to_score2.keys())
		#print comb_id, drug1, drug2, cell_line, len(genes_common)
		values = [ len(genes_common) ]
		if use_expression == "pathway-genes":
		    # Mutual coverage of pathways
		    values_inner = [ 0 ] * len(genes_pathway)
		    for k, gene in enumerate(genes_pathway):
			val = 0
			if gene in genes_to_score1: 
			    val += genes_to_score1[gene]
			if gene in genes_to_score2:
			    val += genes_to_score2[gene]
			values_inner[k] = val
			#if gene in genes_to_score1:
			#    values_inner[k] += 1
			#if gene in genes_to_score2:
			#    values_inner[k] += 1
		    #values.extend(values_inner)
		    values.extend([numpy.mean(values_inner), numpy.max(values_inner)])
		else:
		    #if len(genes_common) == 0:
		    #	print "No genes in common:", comb_id
		    #	continue
		    values_inner = []
		    if use_expression == "no": 
			if comb_to_score is None or len(genes_common) == 0:
			    values_inner = [0]
			else:
			    for gene in genes_common:
				val = comb_to_score[gene] - (genes_to_score1[gene] + genes_to_score2[gene]) / 2
				values_inner.append(val)
			values.extend([numpy.median(values_inner), numpy.max(values_inner)])
		    elif use_expression == "filter":
			genes_common_new = set()
			for gene in genes_common:
			    if gene not in gene_to_idx or cell_line not in cell_line_to_idx:
				continue
			    exp = gexp_norm[gene_to_idx[gene], cell_line_to_idx[cell_line]]
			    if exp < 1:
				continue
			    genes_common_new.add(gene)
			genes_common = genes_common_new
			if comb_to_score is None or len(genes_common) == 0:
			    values_inner = [0]
			else:
			    for gene in genes_common:
				val = comb_to_score[gene] - (drug1_to_score[gene] + drug2_to_score[gene]) / 2
				values_inner.append(val)
			values.extend([numpy.median(values_inner), numpy.max(values_inner)])
		    elif use_expression == "exact":
			for gene in genes_common:
			    if gene not in gene_to_idx or cell_line not in cell_line_to_idx:
				continue
			    val = gexp_norm[gene_to_idx[gene], cell_line_to_idx[cell_line]]
			    values_inner.append(val)
			if len(genes_common) == 0:
			    values_inner = [0]
			values.extend([numpy.median(values_inner), numpy.max(values_inner)])
		    else:
		        raise ValueError("Uknown parameter: %s" % use_expression)
		cell_line_to_values[cell_line] = values
	    combination_to_values[comb_id] = cell_line_to_values
    #cPickle.dump(combination_to_values, open(dump_file, 'w')) 
    return combination_to_values 


def get_top_genes_and_scores(score_file):
    node_to_score = guild_utilities.get_node_to_score(score_file)
    #print node_to_score.items()[:3]
    values = node_to_score.items()
    values.sort(key=lambda x: x[1])
    return node_to_score, values[-int(CONFIG.get("n_top_in_guild")):]


def guild_combinations(drug_to_values, nodes, gexp_norm, gene_to_idx, cell_line_to_idx):
    drugs = drug_to_values.keys()
    flag = False # True #
    z_value = float(CONFIG.get("z_value"))
    for i, drug1 in enumerate(drugs):
	#if drug1 == "NAE":
	#    flag = True
	#if drug1 == "Docetaxel":
	#    flag = False
	#if not flag:
	#    continue
	for j, drug2 in enumerate(drugs):
	    if i < j:
		targets1, smiles1 = drug_to_values[drug1]
		targets2, smiles2 = drug_to_values[drug2]
		targets = targets1 | targets2
		if drug1.lower() < drug2.lower():
		    combination = "%s.%s" % (drug1, drug2)
		else:
		    combination = "%s.%s" % (drug2, drug1)
		combination = "%s.%s" % (drug1, drug2)
		#print combination, targets
		#if len(targets) == max(len(targets1), len(targets2)):
		#    print "Combination redundant!"
		#    continue
		#if len(targets & nodes) == 0:
		#   print "Not in network!"
		#   continue
		if gexp_norm is not None:
		    for cell_line, idx in cell_line_to_idx.iteritems():
			target_to_score = {}
			for target in targets: 
			    if target not in gene_to_idx:
				target_to_score[target] = z_value
			    else:
				val = gexp_norm[gene_to_idx[target], idx]
				target_to_score[target] = z_value if val < z_value else val
			run_guild(combination + ".%s" % cell_line, target_to_score, nodes, qname="all.q")
		else:
		    target_to_score = {}
		    for target in targets: 
			target_to_score[target] = 1
		    run_guild(combination, target_to_score, nodes, qname="all.q")
    return


def guild_drugs(drug_to_values, nodes, gexp_norm, gene_to_idx, cell_line_to_idx):
    z_value = float(CONFIG.get("z_value"))
    for drug, values in drug_to_values.iteritems():
	targets = values[0]
	#print drug, targets
	#if len(targets & nodes) == 0:
	#   print "Not in network!"
	#   continue
	if gexp_norm is not None:
	    for cell_line, idx in cell_line_to_idx.iteritems():
		if cell_line == "647-V":
		    print cell_line, targets
		target_to_score = {}
		for target in targets: 
		    if target not in gene_to_idx:
			target_to_score[target] = z_value
		    else:
			val = gexp_norm[gene_to_idx[target], idx]
			target_to_score[target] = z_value if val < z_value else val
		run_guild(drug + ".%s" % cell_line, target_to_score, nodes, qname="all.q") 
	else:
	    target_to_score = {}
	    for target in targets: 
		target_to_score[target] = 1
	    run_guild(drug, target_to_score, nodes, qname="all.q") 
    return


def run_guild(drug, target_to_score, nodes, qname=None): #targets
    # Create node file
    create_node_file(drug, target_to_score, nodes)
    # Get parameters
    executable_path = CONFIG.get("guild_path") 
    output_dir = CONFIG.get("output_dir") + "/"
    network_lcc_file = CONFIG.get("network_file") + ".lcc"
    output_file = "%s%s.ns" % (output_dir, drug)
    node_file = "%s%s.node" % (output_dir, drug)
    n_repetition = int(CONFIG.get("n_repetition"))
    n_iteration = int(CONFIG.get("n_iteration"))
    # Get and run the GUILD command
    #print strftime("%H:%M:%S - %d %b %Y") #, score_command
    score_command = ' -s s -n "%s" -e "%s" -o "%s" -r %d -i %d' % (node_file, network_lcc_file, output_file, n_repetition, n_iteration)
    if qname is None:
	score_command = executable_path + score_command
	os.system(score_command)
    else:
	#os.system("qsub -cwd -o out -e err -q %s -N %s -b y %s" % (qname, scoring_type, score_command))
	#print "qsub -cwd -o out -e err -q %s -N guild_%s -b y %s" % (qname, drug, score_command)
	print "%s" % (score_command.replace('"', ''))
    return


if __name__ == "__main__":
    main()
    
