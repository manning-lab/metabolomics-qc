# Metabolomics QC for DoD project

Note that this pipeline is for reference rather than for use. Some code will not work out-of-the-box, as it includes code for wrangling data that was in a non-standard format.

The code that is the most reusable is `02_qc_mets.R`, which accepts metabolite measurement data with metabolites as rows, samples as columns.

Data was not pushed to GitHub as it is restricted-access, but examples of what it looked like are given below.

# Scripts
+ `00_install_pkgs.R`: installs required R packages.
+ `01_preprocess_mets.R`: Takes the raw data that was in non-standard format, extracts metabolite metadata, sample metadata, and metabolite measurements. Example below:
  - Raw data:
```
 |         |          |          |      |                           | Date extracted        | MM/DD/YY | MM/DD/YY
 |         |          |          |      |                           | Date injected         | MM/DD/YY | MM/DD/YY
 |         |          |          |      |                           | Column                | 1        | 1
 |         |          |          |      |                           | Analysis order        | 2        | 3
 | Method  | Compound | m/z      | RT   | HMDB ID (*representative) | Metabolite            | TOM###1  | TOM###2
 | HIL-pos | QI221    | 126.1369 | 7.45 | Internal Standard         | valine-d8             | 41368    | 42225
 | HIL-pos | QI197    | 174.1371 | 6.75 | Internal Standard         | phenylalanine-d8      | 136352   | 139558
 | HIL-pos | QI2      | 197.0677 | 3.58 | HMDB11103                 | 1.7-dimethyluric acid | 5141     | 648
```
  - Extracted metabolite metadata:
```
| Method  | Compound | m/z      | RT   | HMDB ID
| HIL-pos | QI221    | 126.1369 | 7.45 | Internal Standard         
| HIL-pos | QI197    | 174.1371 | 6.75 | Internal Standard         
| HIL-pos | QI2      | 197.0677 | 3.58 | HMDB11103                 
```
  - Extracted sample metadata:
```
extraction_date | injection_date | column | analysis_order | sample_id
MM/DD/YY        | MM/DD/YY       | 1      | 2              | TOM###1
MM/DD/YY        | MM/DD/YY       | 1      | 3              | TOM###2
```
  - Extracted measurement data:
```
           TOM###1 | TOM###2
QI221_hn | 41368   | 42225
QI197_hn | 136352  | 139558
QI2_hn   | 5141    | 649
```
+ `02_qc_mets.R`: Performs QC on the formatted measurement data. Steps:
  - Remove metabolite with 0 variance
  - Remove metabolites with >25% missingness
  - Remove samples with >25% missingness
  - Apply various transformations to the data
  - Winsorize to the 5th standard deviation
  - Batch adjustment by median-scaling to a reference batch
+ `03_merge_mets.R`: There were separate files for each LC-MS method (C8-positive, HILIC-positive, Amines, [see here](https://topmed.nhlbi.nih.gov/sites/default/files/TOPMed_CORE_Year3_Broad-BIDMC_metabolomics_methods.pdf)). This script merges the methods together. 

