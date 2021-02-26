python create_histogram.py -gene_list /shared/jdafflon/gene_list_cnn.csv -connectivity_list /shared/jdafflon/thr_networks/corr_coeff_thr_CREBBP.csv -data /shared/jdafflon/cleaned_cnn_data.csv -save_dir NEPDF_data
python train_model.py -save_dir /projects/07e6ba69-8021-439f-901e-4d6226e5b28a/CNNC
