# HDRUK Spring Project: CRUK (CNN gene expression data)

**The paths in the code files refer to the development on the turing DSG system,
you will need to amend them in each case.**

## Pre-processing

The initial step required to run the CNN is to normalise and filter date.

The existing code for pre-processing of the data is in `src/processing`.

This step is performed using `R` and requires the packages in
`src/processing/0_lib_installation.R`. Here is a list.

- `edgeR`
- `DESeq2`
- `tidyverse`
- `pheatmap`
- `EnhancedVolcano`
- `UpSetR`

Here are the main files to look at:

- **`create-norm-data.R`**: normlaise the data it outputs normalised count
  matrix and a log transformed normalised matrix.
- **`DiffExpression_DESeq`**: calculate differential expression either overall
  or per sample.

Some of the other files allow you to calculate log fold change and other
properties. Other files create plots etc.

## Correlation

The output from the correlation network is currently used as the input to the
CNN. The code to get those files are found in the `src/correlation` folder.

# CNN

The CNN code can be found in `src/cnn-histogram`. Two CNN methods were attempted
hence the addition of histogram.

The first step is described in `create_histogram.py`. The rest should be self
explanatory.

Further description and etails on the code can be found in
`src/cnn-histogram/README.md`.