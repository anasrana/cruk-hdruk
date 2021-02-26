# CNN using the histograms

This folder contains all the information necessary to run the CNN model
introduced by [Yuan et. al, 2019](https://www.pnas.org/content/116/52/27151) and
using the data provided by the DSG  challenge.

## Overview of the files:
- `commands_list.sh`: Overview on how to call the main scripts used for the
  analysis with the appropriate arguments. 
- `create_histogram.py`: This scripts takes the log2 normalised data for the
    genes of interest and save the 2D histograms for all gene pais. Note:
    because of memory issues the histograms are saved into 4 different numpy
    files.
    This file is an adaption from the original code from the [CNNC](https://github.com/xiaoyeye/CNNC/tree/master/data)
    package.
- `Jessica_vis_df.csv`: File with the updated network after the training from
  the CNN. Note that not all pair of gene interaction was updated. Only the
  genes that were used for testing were updated.
- `Jessica_visulisation.py`: Script to create the top figure from
  Fig:cnn_results from the report. Note that only 20 randomly selected  genes are
  used for visualisation purposes.
- `Explore histograms.ipynb`: Visualise some of the 2d histograms that were
  created and plot the final heatmap to visualise the results from the CNN.
- `requirements.txt`: List of python packages necessary to run the code
- `train_model.py`: This script is used to train the CNN. 
   Script adapted from the original code from the [CNNC](https://github.com/xiaoyeye/CNNC/tree/master/data)
    package.

