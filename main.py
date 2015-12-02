
from toolbox import configuration, guild_utilities, wrappers
from time import strftime
import os, csv, numpy
try:
    from indigo import Indigo
except:
    print "Indigo not found, drug chemical similarity will not be available!" #!

CONFIG = configuration.Configuration() 


def main():
    # Get network
    network = wrappers.get_network(network_file = CONFIG.get("network_file"), only_lcc = True)
    nodes = set(network.nodes())
    #create_edge_file(network)
    # Get drug info
    drug_to_targets = get_drug_info(network_nodes=nodes)
    #print drug_to_targets.keys()
    # Check individual drugs
    #guild_drugs(drug_to_targets, nodes)
    # Check pairwise drug combinations
    #guild_combinations(drug_to_targets, nodes)
    # Get synergy info
    combination_to_values = get_synergy_info()
    # Get gexp info
    #! gexp_norm, gene_to_idx, cell_line_to_idx = wrappers.get_expression_info(gexp_file = CONFIG.get("gexp_file"), process=set(["z", "abs"]))
    # Check synergy between known pairs
    #check_synergy(combination_to_values, gexp_norm, gene_to_idx, cell_line_to_idx)
    # Get mutation info
    gene_to_cell_line_to_mutation = get_mutation_info(CONFIG.get("mutation_file")) 
    print gene_to_cell_line_to_mutation.items()[:3]
    # Get CNV info
    gene_to_cell_line_to_cnv = get_cnv_info(CONFIG.get("cnv_file")) 
    print gene_to_cell_line_to_cnv.items()[:3]
    # Get pathway info
    #pathway_to_geneids = get_pathway_info(None)
    # Create data file
    return


def get_pathway_info(nodes):
    pathway_to_geneids = wrappers.get_pathway_info(pathway_file = CONFIG.get("pathway_file"), prefix = CONFIG.get("pathway_source"), nodes = nodes)
    print len(pathway_to_geneids) 
    #values = map(lambda x: (x[0], len(x[1])), pathway_to_geneids.items())
    #values.sort(key=lambda x: x[1])
    #print values[-10:]
    genes = pathway_to_geneids["kegg_pathways_in_cancer"]
    genes_merged = set()
    for key, values in pathway_to_geneids.iteritems():
	if key.find("cancer") != -1 and key != "kegg_pathways_in_cancer":
	    print key, len(values)
	    genes_merged |= pathway_to_geneids[key]
	elif key.find("growth") != -1 or key.find("apoptosis") != -1:
	    print key, len(values)
    print len(genes), len(genes_merged), len(genes & genes_merged)
    return pathway_to_geneids 


def check_synergy(combination_to_values, gexp_norm, gene_to_idx, cell_line_to_idx):
    output_dir = CONFIG.get("output_dir") + "/"
    out_file = CONFIG.get("output_file")
    use_expression = int(CONFIG.get("use_expression"))
    f = open(out_file, 'w')
    f.write("comb.id cell.line max.a max.b med mean sd max min syn\n")
    for comb_id, cell_line_to_values in combination_to_values.iteritems():
	drug1, drug2 = comb_id.split(".")
	for cell_line, vals in cell_line_to_values.iteritems():
	    #print comb_id, drug1, drug2, cell_line
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
	    print comb_id, len(genes_common)
	    if len(genes_common) == 0:
		continue
	    comb_to_score = guild_utilities.get_node_to_score(comb_file)
	    values = []
	    if use_expression == 1:
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
	    elif use_expression == 0:
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
	    max_a, max_b, synergy = vals
	    f.write("%s %s %f %f %f %f %f %f %f %s\n" % (comb_id, cell_line, max_a, max_b, numpy.median(values), numpy.mean(values), numpy.std(values), numpy.max(values), numpy.min(values), synergy))
    return


def get_top_genes_and_scores(score_file):
    node_to_score = guild_utilities.get_node_to_score(score_file)
    #print node_to_score.items()[:3]
    values = node_to_score.items()
    values.sort(key=lambda x: x[1])
    return node_to_score, values[-500:]


def guild_combinations(drug_to_targets, nodes):
    drugs = drug_to_targets.keys()
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
		targets1 = drug_to_targets[drug1]
		targets2 = drug_to_targets[drug2]
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


def guild_drugs(drug_to_targets, nodes):
    for drug, targets in drug_to_targets.iteritems():
	print drug, targets
	#if len(targets & nodes) == 0:
	#   print "Not in network!"
	#   continue
	run_guild(drug, targets, qname="all.q")
    return


def get_mutation_info(file_name):
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
    return gene_to_cell_line_to_value


def get_cnv_info(file_name):
    gene_to_cell_line_to_value = {}
    f = open(file_name)
    reader = csv.reader(f, delimiter=",", quotechar='"')
    header = reader.next()
    print header
    header_to_idx = dict((word, i) for i, word in enumerate(header))
    #print header_to_idx
    for row in reader:
	try:
	    ch = row[header_to_idx['chr_GRCh37']]
	except:
	    print row #!
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
    return gene_to_cell_line_to_value


def get_drug_info(network_nodes=None):
    drug_file = CONFIG.get("drug_file")
    drug_to_values = {}
    f = open(drug_file)
    reader = csv.reader(f, delimiter=',', quotechar='"')
    for row in reader:
	drug = row[0]
	targets = set(row[1].split(","))
	formula = row[2]
	if network_nodes is not None:
	    if len(targets & network_nodes) == 0:
		continue
	drug_to_values[drug] = (targets, formula)
	#print drug, targets
    f.close()
    return drug_to_values


def get_synergy_info():
    combination_file = CONFIG.get("combination_file")
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


def get_drug_similarity(drug_to_values, drug1, drug2, method):
    targets1, smiles1 = group_to_values[drug1]
    targets2, smiles2 = group_to_values[drug2]
    targets_common = targets1 & targets2
    d = None
    if method == "target":
	d = len(targets_common) / float(len(targets1|targets2)) # float(min(len(targets1), len(targets2)))
    elif method == "chemical":
	indigo = Indigo()
	m = indigo.loadMolecule(smiles1) 
	m.aromatize()
	fp = m.fingerprint("sim") # sub
	m2 = indigo.loadMolecule(smiles2) 
	m2.aromatize() # Aromatize molecules in case they are not in aromatic form
	fp2 = m2.fingerprint("sim") # Calculate similarity between "similarity" fingerprints
	d = indigo.similarity(fp, fp2, "tanimoto") # tversky
	#print group, group2, "Tanimoto: %s" % (d) 
    else:
	raise ValueError("Uknown method: %s" % method)
    return d


if __name__ == "__main__":
    main()
    
