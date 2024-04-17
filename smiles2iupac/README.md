# smiles_to_iupac.py

This script, "smiles_to_iupac.py", facilitates the conversion of chemical compounds represented in Simplified Molecular Input Line Entry System (SMILES) notation to International Union of Pure and Applied Chemistry (IUPAC) names.

## Usage

**Parameters:**

- `-f FILE`: Specifies the input CSV file containing chemical data.
- `--smiles_col COLUMN_NAME`: Indicates the name of the column in the CSV file where SMILES strings are stored. For example, if SMILES strings are in a column labeled "smiles", specify `--smiles_col "smiles"`.
- `--index_col COLUMN_NAME`: Specifies the name of the column in the CSV file where index numbers are located. For example, if index numbers are in a column labeled "num", specify `--index_col "num"`.

**Example Usage:**

```
python smiles_to_iupac.py -f BBBP5.csv --smiles_col "smiles" --index_col "num"
```

**Example Output:**

| num  | SMILES                                                       | IUPAC_name                                                   | SMILES_reverse                                               | are_same_molecule | tanimoto_similarity | time_elasped       |
| ---- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ----------------- | ------------------- | ------------------ |
| 1    | [Cl].CC(C)NCC(O)COc1cccc2ccccc12                             | 1-naphthalen-1-yloxy-3-(propan-2-ylamino)propan-2-ol;chloride | CC(C)NCC(COC1=CC=CC2=C1C=CC=C2)O.[Cl-]                       | False             | 0.9642857142857143  | 10.536948680877686 |
| 2    | C(=O)(OC(C)(C)C)CCCc1ccc(cc1)N(CCCl)CCCl                     | tert-butyl4-[4-[bis(2-chloroethyl)amino]phenyl]butanoate     | CC(C)(C)OC(=O)CCCC1=CC=C(C=C1)N(CCCl)CCCl                    | True              | 1.0                 | 1.2848384380340576 |
| 3    | c12c3c(N4CCN(C)CC4)c(F)cc1c(c(C(O)=O)cn2C(C)CO3)=O           | 7-fluoro-2-methyl-6-(4-methylpiperazin-1-yl)-10-oxo-4-oxa-1-azatricyclo[7.3.1.05,13]trideca-5(13),6,8,11-tetraene-11-carboxylicacid | CC1COC2=C(C(=CC3=C2N1C=C(C3=O)C(=O)O)F)N4CCN(C)CC4           | True              | 1.0                 | 2.919210433959961  |
| 4    | C1CCN(CC1)Cc1cccc(c1)OCCCNC(=O)C                             | N-[3-[3-(piperidin-1-ylmethyl)phenoxy]propyl]acetamide       | CC(=O)NCCCOC1=CC(=CC=C1)CN2CCCCC2                            | True              | 1.0                 | 1.497406244277954  |
| 5    | Cc1onc(c2ccccc2Cl)c1C(=O)N[C@H]3[C@H]4SC(C)(C)[C@@H](N4C3=O)C(O)=O | (2S,5R,6R)-6-[[3-(2-chlorophenyl)-5-methyl-1,2-oxazole-4-carbonyl]amino]-3,3-dimethyl-7-oxo-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylicacid | CC1=C(C(=NO1)C2=C(C=CC=C2)Cl)C(=O)N[C@@H]3C(=O)N4[C@@H](C(=O)O)C(C)(C)S[C@H]34 | True              | 1.0                 | 3.3189845085144043 |





# Install 

```
conda create --name STOUT python=3.8 
conda activate STOUT
pip install STOUT-pypi
```





# References

1. [Kohulan/Smiles-TO-iUpac-Translator](https://github.com/Kohulan/Smiles-TO-iUpac-Translator)

