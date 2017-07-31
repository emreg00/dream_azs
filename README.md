
# Understanding molecular bases of drug response and drug synergy
Understanding synergistic effects of drugs is key to develop effective intervention strategies targeting diseases (such as AD or T2D or both) and provides unprecedented opportunities to repurpose existing drugs. The AstraZeneca-Sanger Drug Combination Prediction DREAM Challenge provides a rich data source aiming to understand the synergistic drug behavior based on pretreatment data and spans cell viability data over 118 drugs and 85 cancer cell lines (primarily colon, lung, and breast). In collaborating with Dr. Baldo Oliva's group at GRIB, UPF-IMIM, we have been working on identifying the effects of confounding factors in the data set such as dosage and genetic background of the cell lines and developing algorithms that can predict the individual and synergistic effects of drugs. 

The first challenge has two subtasks: predicting drug synergy *(i)* using mono synergy data *(ii)* without using mono synergy and gene expression data. The participants are free to use any other data source (such as cell line data, gene mutation, drug target information provided in the challenge or external data sets) and submit their predictions in 3-4 rounds, which is followed by a final round. The second challenge requires making predictions for drug combinations and cell lines for which no previous training data available (making it hard to build a machine learning predictor). 

Check challenge info and timelines at https://www.synapse.org/#!Synapse:syn4231880/wiki/235652

Next deadline (for both tasks): March, 14th 2016 (Final round)

<For the first round of the challenge, we have build machine learning models to predict the synergy of drugs for both of these tasks and choice the best performing models to submit predictions. Among various machine learning models, we found a combination of bootstrapped and ensemble tree-based predictors achieved best performance on the training data set for. 

To improve the prediction performance we have incorporated mutation data (of drug targets in a given cell line) and interactome based contribution of the drug combination compared to the effect of drugs separately. To assess interactome based contribution of a drug or combination (characterized by a set of targets), we have used GUILD, a network-based functional prioritization tool. 

Interestingly, using GUILD, only the predictions for subtask *(ii)* improved but not for subtask *(i)*. We suspect this is due to the mono therapy response data describing the synergy best and addition of new features (such as the ones based on expression, mutation, interactome) potentially causing the predictor to overfit to the training data set.>

## Data overview

Challenge training data consists of 2199 samples providing information on 169 drug pairs over 85 cell lines.


```R
source("dream.R")
parameters = initialize()
dat = overview(parameters, summarize=F)
```

    [1] "Number of samples: 2199"
    [1] "Number of drug pairs: 169, number of cell lines: 85"
    [1] "Features: CELL_LINE, COMPOUND_A, COMPOUND_B, MAX_CONC_A, MAX_CONC_B, IC50_A, H_A, Einf_A, IC50_B, H_B, Einf_B, SYNERGY_SCORE, QA, COMBINATION_ID"


## Data cleaning, preprocessing and imputation

- Filter samples with low quality (404 samples):
$$QA < 1$$


- Filter samples with low sensitivity (3 cell lines) based on the observation that higher Einf correlates with lower synergy. Define Einf of the drug pair A,B as follows:
$$min((Einf_A + Einf_B) / 2)$$ 


- Filter correlated features (None):
$$ PCC > 0.75 $$

Both min (-588.221) and max (6737.175) synergy instances have low quality. After filtering synergy scores range between -179 and 237.

- Filter instances in which CNV values are NA

- Scale all the features using z-score transformation (centered and scaled by standard deviation)


- Impute missing values using k-nearest-neighbor ($k = 5$)


- (For challange 2) The synergy values are categorized as follows:
$$
category = \{
\begin{array}{cl}
1, & if~synergy > 30 \\\
0, & otherwise
\end{array}
$$
The negative cases ($P = \{x: x=0\}$, where x corresponds to the synergy category of the instance) are more abundant than positive cases ($ N = \{x: x=1\}$, ~80% vs ~20%). Accordingly, we balanced the data set such that $|P| = |S|$.


```R
dat = filter(dat, cutoff=40)
```

    [1] "Correlated features:"
    [1] "Correlation between einf.min and syn.med: -0.235249 0.031231"
    [1] "Insensitive cell lines: 22RV1, KU-19-19, VCaP"
    [1] "Number of samples with QA < 1: 404, Einf > 40: 7"


## Target prediction

We use SMILES and target information of all drugs in DrugBank to predict targets of the drugs in the data set (Drug_info_release). The Tanimoto coefficient between aromitized SMILES of drugs is used to define similarity (Tanimoto coefficient cutoff=0.6). We calculate a Fisher's based enrichment score to find the overrepresented targets between similar SMILES (FDR cutoff=20% after Benjamini Hochberg multiple hypothesis test correction procedure). We also manually checked PubChem for targets and used HitPick server to find drug targets. We considered SEA for target prediction but then did not use it since (i) it did not offer a batch run mode (for multiple querries at a time) and (ii) its predictions seemed overlapping (~promiscous). Predicted targets are used as seeds while running GUILD, accordingly only taken into account in GUILD-based features. The drug targets span 171 genes (after manual curation) and 212 genes after SMILES based prediction.


## Feature definition

***Baseline prediction***

- Monotherapy response based
    * max concentration
    * viability at max kill
    * IC50 
    * slope of the fit to the dose response curve
    
***Expression based***

Expression values of each gene in a given cell line are converted to z-scores using the average/s.d. of the distribution of the gene's expression over all the cell lines.

- The average gene expression of the targets of two drugs in the cell line
    * $gexp = abs(gexpA) + abs(gexpB)$, where for each cell line where gexpA is the average expression of A's targets T and
    
$$ gexp(T, cell) = median \{ t \in targets(T) | E(t, cell) \} $$ where $E$ is the gene expression matrix, $T$ is the drug tested in a given $cell$ line.

- For each drug target, the sum of expression in the cell line if it is the target of one of the two drugs in the combination (171 features)
    * $gexpT = \sum_{T \in all~targets} X_T * gexp(T, cell)$, where X_T is the number of drugs in the combination for which T is a target. "All targets" are all 212 the targets in the data set (Drug info release + SMILES based target prediction).

***Mutation based***

- The average mutation score of the targets of the two drugs in the cell line
    * $mut = abs(mutA - mutB)$ , where for each cell line where mutT is
    
$$ mut(T, cell) = median \{ t \in targets(T) | M(t, cell) \} $$ where $M$ is the mutation, $T$ is the drug tested in a given $cell$ line. Genes are assigned mutation score based on the "Description" field in the annotation file (0 if the mutation is silent or of unknown impact, 2 if the mutation is associated to cancer with respect to FATHMM prediction and 1 otherwise). Impute missing values using k-nearest-neighbor ($k = 5$).

- For each drug target, the sum of mutation scores in the cell line if it is the target of one of the two drugs in the combination (171 features)

***Copy numbre variation based***
- The difference between the copy number variation (max in case of difference in diploid copy) of the targets of the two drugs in the cell line
    * $cnv = abs(cnvA - cnvB)$ , where for each cell line where cnvT is
    
$$ cnv(T, cell) = \{ t \in targets(T) | C(t, cell) \} $$ where $C$ is the copy number, $T$ is the drug tested in a given $cell$ line. Impute missing values using k-nearest-neighbor ($k = 5$).

- For each drug target, the sum of CNV scores in the cell line if it is the target of one of the two drugs in the combination (171 features)

***Interactome based***

- The network-impact score distribution of the genes in the overlap between top 500 genes in GUILD-based prioritization of drug targets of A and B, respectively. 
    * guild.common (number of common genes)
    * guild.med (mean of the distribution of the network impact)
    * guild.max (mean of the distribution)

The network-impact is calculated as
$$ impact(A,B) = GUILD({A,B}) - (GUILD(A) + GUILD(B)) / 2 $$ where $GUILD(X)$, is the GUILD scores of the genes when genes in X are used as seeds. Top scoring 500 genes common in $GUILD(A)$ and $GUILD(B)$ are considered to calculate the impact score distribution.

- Target degree based
    + difference ($abs(kA-kB)$) 
    + max ($max(kA, kB)$) 
    + min ($min(kA, kB)$)

- Distance between targets: d(A, B)
    
***Drug similarity based***

- If the drugs are similar, the effect is expected to be synergistic (i.e. Loewe additivity)
    * sim.target: common targets
        $$ sim(A, B) = \frac{T(A) \cap T(B)}{T(A) \cup T(B)}  $$
        
    * sim.chemical: chemical formula similarity, calculated using Tanimota similarity coefficient (Jaccard index of molecular fingerprints). 

***KEGG pathways***

- Cancer related from  KEGG pathways. These pathways are "pathways in cancer", "aminoacyl-tRNA biosynthesis", "MAPK signaling pathway", "NF-kappa B signaling pathway". For genes in these pathways,
    * involvement of drug targets in these pathways (kegg.in, 2: targets of both drugs in combination are in the pathway, 1: only targets of one are in the pathway, 0: none of the targets are in the pathway)
    * gene expression (kegg.gexp.med and kegg.gexp.max: the median and max of the distribution)
    * mutation (kegg.mut.med and kegg.mut.max)
    * CNV (kegg.cnv.med and kegg.cnv.max)
    
***Cancer genes***

- COSMIC genes from http://cancer.sanger.ac.uk/census/ (572 genes).
    * involvement of drug targets in these pathways (cosmic.in)
    * gene expression (cosmic.gexp.med and cosmic.gexp.max)
    * mutation (cosmic.mut.med and cosmic.mut.max)
    * CNV (cosmic.cnv.med and cosmic.cnv.max)

***Categorized features***

- For each pathway gene (74 genes in COSMIC + KEGG pathways above), we create an individual feature denoting whether the drug target is that pathway gene
- For each drug target (171 genes w/o predictions), we create an individual feature for the following
    * .gexp (expression for that target on that cell line)
    * .mut (0/1/2 for that target on that cell line)
    * .cnv (average CNV for that target on that cell line)
    * .zygosity (from CNV, "H": 0, "LOH": 1, "0": 2 for that target on that cell line)
    
***Combined***
- Features used in the final models (unless otherwise stated below)
    + gexp
    + mut
    + cnv
    + guild.common
    + guild.med
    + guild.max
    + sim.target
    + sim.chemical
    + kegg.in
    + cosmic.in
    + .gexp(171)
    + .mut(171)
    + .cnv(171)
    + .zygosity(171)
    + .pathway(74)
    + k.diff
    + k.min
    + k.max
    + dAB
    
## Feature definition and prediction models

- For ~~challange 1~~ both challenges, the best performing model was using RandomForest and Generalized Boosted Regression Models with the combination of features above

- ~~For challange 2, generalized linear model was used~~



```R
results(parameters)
```


    Error in eval(expr, envir, enclos): could not find function "results"



\# 8 above is under consideration for Round 4

Round 3 is better than 4. Could be due to
- drug target info change
- expression abs / definition change
- cnv definition change
- used combination of features

## General considerations
* The training set performance does not necessarily correlate with evaluation set performance (e.g., [R3] outperforms [R2] in the above table on DREAM's evaluation)
* The predictor is more robust if the rf and gbm have similar training set performance (rather than combined performance). That is the predictor should both have good combined performance and the performance of the classifiers, ideally, should be closer to each other.
* The scoring script provided by organizers is a better descriptive of the performance
* Using full training data set (as opposed to 70% split used for training & validation), slightly improves the performance
* Non-cell-line specific GUILD features outperformed cell line specific features (without considering cell line specific combinations, results currently not generated due to the large number of runs $(120/2) * 119 * 83$ GUILD runs)
* Challange 2 requires making predictions for cell lines and drug combinations for which no (or little) training data exists. This is a modeling task rather than machine learning task (as highlighted by the organizers). Potentially the drug response on a given cell line can be modeled using external data and a predictor of synergy based on the combination of these responses can be built. At the current stage, we treat this as a machine learning problem, building a binary classifier using the features in challange 1.
* Combination of predictors improves the performance substantially

## Final predictors and confidence assignment

>***Predictor for challenge 1 subtask 1*** 
The best performing predictor using the response data and the features above achieves an accuracy (assessed by correlation between predicted and observed synergy scores) of 0.4 on the training set and 0.2 on the evaluation set. Note that the value on the training set fluctuates depending on the folds used in cross validation (+/- 0.1). 

>***Predictor for challenge 1 subtask 2*** 
Achieves an accuracy (assessed by correlation between predicted and observed synergy scores) of 0.4 on the training set and 0.15 on the evaluation set.
    
>***Predictor for challenge 2*** 
This challenge requires predicting drug combinations and cell lines for which no previous training data is available, thus makes it harder to find features that would work over all the test data (due to the missing values). The predictor achieves a F1 score of 0.64 on the training data set and 0.36 on the evaluation set.

>***Confidence assignment***
We observed that the predictions tend to fail for higher synergy scores, accordingly we defined the following confidence scoring:
$$confidence = 1 - abs(synergy) / max(abs(synergy))$$

## DREAM evaluation 

>The global correlation values from DREAM (assessed by real values in the test set) are substantially lower than the correlation values in the training set (assessed by model development on 70% of the training data using cross validation and validation using 30% of the data).
* Challenge 1 subtask 1: 0.18 (global), 0.25 (mean of top 10-20-30%), 
* Challenge 1 subtask 2: 0.15 (global), 0.21 (mean of top 10-20-30%)
* Challenge 2: 0.36 (global F1_20), 0.21 (mean 1-way row ANOVA of top 10-20-30%)
* Overall ranking: among 20% of all submissions for challenge 1, among %70 for challange 2

The table below shows global and top (corresponding to highest confidence 30%) performance for various tasks as well as the values for the best ranking submission. For challange 1 the global and top measures are global correlation and mean correlation over top 10-20-30%, respectively. For challenge 2, the global and top measures are F1 score (using Combofit score > 20) and mean 1-way row-wide ANOVA, respectively.


```R
results(parameters, leaderboard=T)
```




<table>
<thead><tr><th></th><th scope=col>Challange</th><th scope=col>X14A41.global</th><th scope=col>X14A41.top</th><th scope=col>Max.global</th><th scope=col>Max.top</th><th scope=col>Rank..top.based.</th><th scope=col>Submissions</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>1A</td><td>0.18</td><td>0.25</td><td>0.27</td><td>0.37</td><td>21</td><td>139</td></tr>
	<tr><th scope=row>2</th><td>1B</td><td>0.15</td><td>0.21</td><td>0.25</td><td>0.32</td><td>26</td><td>112</td></tr>
	<tr><th scope=row>3</th><td>2</td><td>0.36</td><td>0.04</td><td>0.48</td><td>0.27</td><td>56</td><td>76</td></tr>
</tbody>
</table>




## TODO

- ~~Enrich drug target data with SEA and HitPick predictions, filter no target drugs~~
- ~~Check the labeling issue in challenge 2 (1:0, 2:1)~~
- ~~Consider confidence scoring based on the sensitivity of the cell lines (Einf on the training cell lines or performance on the training data set)~~
- ~~Try combination of all features (bottom of table)~~
- ~~Check scoring evaluation scripts by DREAM~~
- ~~Correct combination name bug~~
- ~~Rerun GUILD~~
- ~~Revise existing featuers (gexp / cnv)~~
- ~~Redefine monotherapy response features (difference of A-B)~~
- ~~Check different models / consider only using RF~~
- Use several features proposed in Sun et al., 2015, Nat Comms "Combining genomic and network characteristics for extended capability in predicting synergistic drugs for cancer" http://www.nature.com/ncomms/2015/150928/ncomms9481/full/ncomms9481.html
    + MI GO BP
    + Unrelated pathway ratio
    + ~~distance between targets in PPI~~
    + ~~degree & centrality in PPI~~
- Add methylation-based features
- GUILD top 500 vs 1000
- Consider response prediction from molecular features to be fed to the synergy predictor
- Incorporate external data / synergy modeling (for challange 2) 



```R

```


```R

```
