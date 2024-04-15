This is the Repository for the planned paper

'Starshaped Subgroup Discovery with Uniform Generalization Guaranteess' 

The implementation of the subgroup discoveries (including statistical tests and regularization) is done with the R-package 'oofos' ( https://github.com/schollmeyer/oofos )

There are different applications from different fields: Social science, psychometrics, biology, machine learning (from which only one will later be contained in the paper):

    Dimensions of social justice ( R Code :...... )
    Item impact in Cognitive Diagnosis Models ( R Code: Application_example_cognitive_diagnosis_modeling.R )
    IRT dataset anxiety ( http://openpsychometrics.org/_rawdata/ )
    Gene expression data set(s)
    UCI dataset credit_g ( from https://www.openml.org/search?type=data&sort=runs&id=31&status=active )


The gene expression datasets can be found under:
https://hemberg-lab.github.io/scRNA.seq.datasets/

The gene filter is described here:

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410170/


Documentation of the files:

additional_functions.R: Code with algorithm for the computation of the empirical Rademcahercomplexity

additional_functions_for_gene_expression_data.R: additional code for the import and preprocessing of gene expression data

application_credit_dataset.R: Code for the starshaped subgroup discovery for the credit dataset

dataanalysis_treutlein.R: Code for the analysis of the gene expression dataset (Treutlein *) Used stylized beweenness is the geometry based stylized betweenness (which uses angles)

dataanalysis_treutlein_absb.R: Code for the analysis of the gene expression dataset (Treutlein *). Used stylized beweenness is the attribute based stylized betweenness (which counts how many attribues violate a formal implication)
relation_VC_dimension_rademacher-complexity.R: Code for the visualization of the relation between the VC dimension and the Rademacher Complexity
used_stylized_betweenness_functions.R: Implementation of the different stylized betweenness relations



Outdated files and folders:

versuch_gs_rademacher_complexity.R
dataanalysis_yan.R
results_credit_data
results_treutlein
results_yan
versuch_gs_rademacher_complexity.R
