library(data.table)

smiles_uniform_1 = readRDS("~/Documents/Epochs/source_inject/smiles_uniform/smiles_uniform_1/smiles_uniform_mock_prospect_samples.RDS")
smiles_uniform_2 = readRDS("~/Documents/Epochs/source_inject/smiles_uniform/smiles_uniform/smiles_uniform_mock_prospect_samples.RDS")

smiles_uniform_1_trim = smiles_uniform_1[lengths(smiles_uniform_1) != 0]
smiles_uniform = c(smiles_uniform_1_trim, smiles_uniform_2)

saveRDS(
  smiles_uniform,
  "~/Documents/Epochs/source_inject/smiles_uniform/smiles_uniform_mock_prospect_samples.RDS"
)