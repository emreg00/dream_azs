# Understanding molecular bases of drug response and synergy: An ensemble learning approach using diverse molecular and interactome-based signatures

##Summary Sentence
To build a predictor quantifying synergy between two compounds, we defined molecular and interactome based features and combined Random Forest and Stochastic Gradient Boosting algorithms.


##Background/Introduction
Understanding synergistic effects of drugs is key to develop effective intervention strategies targeting diseases and provides unprecedented opportunities to repurpose existing drugs. 
Existing studies have suggested defining genomics, transcriptomics, pathway involvement and interactome-based features to explain drug synergy (Sun et al., 2015). 
During the initial rounds of the challenge, we have observed that even a simple machine learning predictor using monotherapy response data performs as good as the average submitte models, suggesting that the key tasks in the challenge are (i) defining descriptive features and (ii) building robust learning models. 
The latter is especially important, given the concerns on the noise and reproducibility of drug synergy data.
Following these ideas, we have defined molecular, pharmacological, pathway-level and interactome-based features to characterize the drug effect. 
In addition to using cell line specific gene expression, methylation, mutation, copy number variation (CNV), zygocity, and compound similarity data, we have incorporated a network-based impact score using a tool we have developed recently (Guney and Oliva, 2012). 
We have also included tissue information of cell lines and pathway involvement of targets of the compounds.
After data cleaning and normalization (i.e. removing insensitive cell lines and missing data points, imputation etc...), we have systematically tested several machine learning approaches and finally built a combined Random Forest and Gradient Boosted Machine predictor.


##Methods

*Gene expression and methylation data*

Methylation ilse M-values are mapped to genes closest to the genomic location (used R TxDB and GRanges). Gene expression and methylation values for each gene are normalized over all cell lines (converted to z-scores using the average/s.d. of the distribution of the gene's expression over all the cell lines).

*Pathway information*

We retrieved cancer related from KEGG. These pathways are "pathways in cancer", "aminoacyl-tRNA biosynthesis", "MAPK signaling pathway", "NF-kappa B signaling pathway",  "cell cycle", "P53 signaling pathway", "apoptosis", "TGF-beta signaling pathway". All the genes in these pathways are simply referred as "kegg". We also used COSMIC genes (572 genes downloaded from http://cancer.sanger.ac.uk/census/).

*Target curation and prediction*

We use SMILES and target information of all drugs in DrugBank to predict targets of the drugs in the data set (Drug_info_release). The Tanimoto coefficient between aromitized SMILES of drugs is used to define similarity (Tanimoto coefficient cutoff=0.6). We calculate a Fisher's based enrichment score to find the overrepresented targets between similar SMILES (FDR cutoff=20% after Benjamini Hochberg multiple hypothesis test correction procedure). We also manually checked PubChem for targets and used HitPick server to find drug targets. We considered SEA for target prediction but then did not use it since (i) it did not offer a batch run mode (for multiple queries at a time) and (ii) its predictions seemed overlapping (~promiscuous). Predicted targets are used as seeds while running GUILD, accordingly only taken into account in GUILD-based features. The drug targets span 171 genes (after manual curation) and 212 genes after SMILES based prediction.

*Network-based impact of drugs*

We used a network-based prioritization tool GUILD, to quantify the gene neighborhood affected in the interactome. We used the known drug targets as seeds for each drug and ran NetScore algorithm with the default parameters (r=3, i=2, seed score=1, non-seed score=0.01). 
We considered top 500 genes in GUILD-based prioritization of drug targets as the "affected neighborhood" and checked the network-impact score distribution of the genes in this neighborhood.
We used the integrated human interactome in a recently published study (Menche et al., 2015).

*Data cleaning, preprocessing and imputation*
- Filter samples with low quality ($QA < 1$, removing 404 samples, all numbers on ch1_training_set, unless otherwise stated).
Both min (-588.221) and max (6737.175) synergy instances have low quality. After filtering synergy scores range between -179 and 237.

- Filter samples with low sensitivity (3 cell lines) based on the observation that higher Einf correlates with lower synergy. Define Einf of the drug pair A,B as 
$Einf = min((Einf_A + Einf_B) / 2)$ and keep cell lines with $Einf \le 40$. 

- Filter correlated features ($PCC > 0.75$).

- Filter instances in which CNV values are NA.

- We expanded the drug information in the release file using DrugBank, PubChem and HitPick resources.

- Impute missing values using k-nearest-neighbor ($k = 5$).

- Scale all the features using z-score transformation (centered and scaled by standard deviation).

- The synergy values are categorized as follows (for challenge 2):
$$
category = \{
\begin{array}{cl}
1, & if~synergy > 20 \\\
0, & otherwise
\end{array}
$$
The negative cases ($P = \{x: x=0\}$, where x corresponds to the synergy category of the instance) are more abundant than positive cases ($ N = \{x: x=1\}$, ~80% vs ~20%). Accordingly, we balanced the data set such that $|P| = |S|$.


*Feature definition*

- 4 features based on the monotherapy response of two drugs A and B. For each of the following feature, a single value is obtained using $value(A,B) = abs(value(A) - value(B))$:
    * conc, max concentration 
    * einf, viability at max concentration
    * ic50, concentration at 50\% viability 
    * h, slope of the fit to the dose response curve
    
- 2 features based on the expression values of each gene in a given cell line:
    * gexp, the average gene expression impact of the targets of two drugs in the cell line ($gexp = gexpA.amed * gexpB.amed$, where $gexpX.med$ is the median expression of $X$'s targets in the cell line $X$ is tested: $ gexp.amed(X, cell) = median \{ T \in targets(X) | abs(E(T, cell)) \} $ where $E$ is the gene expression matrix, $X$ is the drug tested in a given $cell$ line)
    * gexp.diff, the gene expression difference between targets of A and B ($gexp = abs(gexpA.amed * gexpB.amed)$)

- 171 features (.g*) based on the expression values of each drug target. The sum of expression in the cell line if it is the target of one of the two drugs in the combination ($gexpX = \sum_{T \in all~targets} N_T * gexp(T, cell)$, where $N_T$ is the number of drugs in the combination for which $T$ is a target. "All targets" are all the targets in the data set (Drug info release)

- 2 features based on the methylation values of each gene in a given cell line, similar to gene expression based values above
    * met, the average gene expression impact of the targets of two drugs in the cell line 
    * met.diff, the methylation difference between targets of A and B 

- 171 features (.e*) based on the methylation values of each drug target similar to gene expression based values above

- 1 feature (mut) based on the average mutation score of the targets of the two drugs in the cell line ($mut = abs(mutA - mutB)$ , where for each cell line where mutT is
$ mut(X, cell) = median \{ T \in targets(X) | M(T, cell) \} $ where $M$ is the mutation, $X$ is the drug tested in a given $cell$ line. 
Genes are assigned mutation score based on the "Description" field in the annotation file (0 if the mutation is silent or of unknown impact, 2 if the mutation is associated to cancer with respect to FATHMM prediction and 1 otherwise). 

- 171 features (.m*) based on mutation profiles of each drug target. The sum of mutation scores in the cell line if it is the target of one of the two drugs in the combination

- 1 feature based on CNV (cnv). The difference between the copy number variation (max in case of difference in diploid copy) of the targets of the two drugs in the cell line ($cnv = abs(cnvA - cnvB)$, where for each cell line cnvX is $ cnv(X, cell) = \{ T \in targets(X) | C(T, cell) \} $$ where $C$ is the copy number, $X$ is the drug tested in a given $cell$ line. 

- 171 features (.c*) for each drug target; the sum of CNV scores in the cell line if it is the target of one of the two drugs in the combination 

- 171 features (.z*) for each drug target; the zygocity of the drug targets ("H": 0, "LOH": 1, "0": 2 for that target on that cell line)

- 3 features based on the impact on the interactome.
The network-impact is calculated as $ impact(A,B) = GUILD({A,B}) - (GUILD(A) + GUILD(B)) / 2 $ where $GUILD(X)$, is the GUILD scores of the genes when genes in X are used as seeds. Top scoring 500 genes common in $GUILD(A)$ and $GUILD(B)$ are considered to calculate the impact score distribution.
    * guild.common, number of common genes among top 500 genes
    * guild.med, median of the distribution of the network impact
    * guild.max, maximum of the distribution

- 4 features based on the targets centrality and closeness in the interactome:
    + k.diff, degree difference ($abs(kA-kB)$) 
    + k.max, ($max(kA, kB)$) 
    + k.min, ($min(kA, kB)$)
    + d, average shortest path distance between targets ($d(A, B)$)
    
- 2 features based on drug similarity.If the drugs are similar, the effect is expected to be synergistic (i.e. Loewe additivity)
    * sim.target, common targets ($ sim(A, B) = \frac{targets(A) \cap targets(B)}{targets(A) \cup targets(B)}  $
    * sim.chemical, chemical formula similarity, calculated using Tanimota similarity coefficient (Jaccard index of molecular fingerprints)

- 2 features based on COSMIC and KEGG pathway involvement 
    * cosmic.in, involvement of drug's targets among COSMIC cancer genes ($inv(A,B) = ||targets(A) \cap cosmic|| + ||targets(B) \cap cosmic||$)
    * kegg.in, involvement of drug's targets among KEGG pathways listed above ($inv(A,B) = ||targets(A) \cap kegg|| + ||targets(B) \cap kegg||$)

- 8 features (.p*) based on the pathway overlap. For 8 pathway from kegg, the sum of the number of drug targets for A and B within that pathway
    
- 2 features based on the median synergy per cell line and drug combination
    * cell.med, median synergy value in the training set per each cell line
    * comb.med, median synergy value in the training set per each drug combination

- 8 features (.t*) based on the tissue, binary (0/1) value corresponding to the tissue the cell line is coming from

*Confidence assignment*
We defined the following confidence scoring: $confidence = 1 - abs(synergy) / max(abs(synergy))$, based on the observation that higher synergy scores tend to be less reliable.

*Prediction model and implementation*
We have used python to parse and process most of the data, created a file containing features and then used R to filter, further process (if need be), build the machine learning model, validate and output the resulting models. 
For data processing and analysis we have used numpy package in python and caret, reshape2, plyr, TxDB and GRanges packages in R. 
For building prediction models we have used rf and gbm packages in R.
A repeated cross-validation approach is used to build the model on the training data set (10-folds, 3 repetitions). 
We used 70% of the data to train the model and the remaining 30% to validate it (assess accuracy using correlation and organizer provided metrics).
For the predictors we used the default parameters in caret (which includes built-in feature selection and parameter optimization).
For replicability we used random seed value 142341 in R.
We used the full training data set in the submission (as opposed to 70% split used for training & validation),

- For challenge 1A, we used all the features above except mut, k.diff, k.max, k.min, d, .m*, .z*
- For challenge 1B, we used all the features in 1A except conc, einf, ic50, h, gexp, met, .g*. .e*
- For challenge 2, we used all the features in 1A but used only Random Forest classifier (instead of RF + GBM)


##Conclusion/Discussion                                    

The first challenge had two subtasks: predicting drug synergy *(i)* using mono synergy and genomics data *(ii)* without using mono synergy, gene expression and methylation data. 
The participants could use any other data source (such as cell line data, gene mutation, drug target information provided in the challenge or external data sets). 
The second challenge requires making predictions for drug combinations and cell lines for which no previous training data available (making it hard to build a machine learning predictor). 
The predictor had near-perfect accuracy on the training data set (assessed via cross-validation) but we suspect that this is mostly due to overfitting. In order to minimize the effect of overfitting, we have submitted a model incorporating less features to challenge 1B and thus, expect it to be more robust.
Based on our experience during the challenge, we summarize the observations we have made below.

*Observations*
- Incorporating more features do not necessarily improve prediction performance, while increasing the risk of overfitting
- Interestingly, mono therapy response data do not necessarily improve the prediction accuracy (which is also reflected by many teams having worse score in 1A compared to 1B)
- Drug similarity, target CNV profile, pathway involvement and network-based features are crucial and robust features contributing to the accuracy of the prediction model
- Combining RF and GBM models tend to result higher accuracy than using them individually


##References
* AZ-Sanger Drug Combination Prediction DREAM Challenge (#!Synapse:syn4231880)
* Sun et al., 2015, Nat Comms, "Combining genomic and network characteristics for extended capability in predicting synergistic drugs for cancer"
* Liu et al., 2013, Bioinformatics, "HitPick: a web server for hit identification and target prediction of chemical screenings"
* Guney E and Oliva B, 2012, PLoS ONE "Exploiting protein-protein interaction networks for genome-wide disease-gene prioritization" 
* Menche et al., 2015, Science, "Uncovering disease-disease relationships through the incomplete interactome"
* [GUILD standalone tool](sbi.imim.es/GUILD.php)
* [Caret package](http://topepo.github.io/caret/index.html)
* Modified drug target information (curated and predicted targets): [Drug_info_release.csv.predicted](#!Synapse:syn5757879)
* R and Python code used in this project: [14A41.R](#!Synapse:syn5757878) and [14A41.py](#!Synapse:5757877)


