# {moduleScoreR}

## Description
This is a simple package that, given an object with RNAseq values (which can be
in a data.frame, DESeqDataObject from {DESeq2}, or a DGEList from {edgeR}) and a
list of gene modules, will assign a score based on the eigenvalues of the first
principle component of the genes for each module.  By default, {moduleScoreR}
uses the `rsvd` function from the {rsvd} package to calculate PC1, but `svd`
from {base} R and `nmf` from the {NMF} package can also be used.

Gene modules need to be in the form of a named list of lists, i.e.
```
$MODULE_A
 [1] "GENE_1" "GENE_2" "GENE_3"

$MODULE_B
 [1] "GENE_4" "GENE_5" "GENE_6"
```
Alternatively, one can prepare the necessary list using `prepGMT` to format a
gmt file from the [Molecular Signatures Database](https://www.gsea-msigdb.org/gsea/msigdb/)
or from a two-column data.frame, with the module names in a column named `term`
and the genes for each module in a column named `gene`

## Installation
```
devtools::install_github("milescsmith/moduleScoreR")
```

## Usage

Basic usage takes the form of:
```
module_list <- prepGMT("path/to/gmt_file")
scoreEigengenes(
  object      = obj,
  module_list = module_list,
  score_func  = "rsvd"
)
```
Optional arguments include:

  - `normalize_func`: choose between using `vst` and `rlog` for DESeqDataSets
    and `cpm` and `rpkm` for DGELists (note, though, that `rpkm` requires gene
    lengths to be present in the DGEList) for preparing expression levels
  - `return_self`: if TRUE, a copy of the object containing expression values
    will be returned with the scores present in the metadata component of the
    respective objects
