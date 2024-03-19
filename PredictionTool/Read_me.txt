Binding energy predictions for MHCI using linear and global epistasis model.

The folder 'PredictionTool' contains energy matrices, coefficients and a code to compute binding energies for the desired alleles and peptides.


## PREDICTION TOOL CONFIGURATION
1. Within the 'Prediction_MHC.py' file, select the desired allele among those available: '0201' for HLA-A*0201, '1101' for HLA-A*1101, '0702' for HLA-B*0702 and 'DRB10101' for HLA-DRB1*0101.

2. Import peptides, better as pandas DataFrames. An example file 'test_out.csv' is already provided for comparison. The only accepted peptide length is 9 aa. For MHCII, this would give and estimate of the binding core energy.

3. Choose model, if linear or quadratic

Individual energy matrices are provided in the folder 'Tables'. 

