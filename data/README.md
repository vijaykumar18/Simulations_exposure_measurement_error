# Data

This folder contains the preterm birth (PTB) data used in the simulation study.

## File

| File | Description |
|------|-------------|
| `zipall_1719.sas7bdat` | ZIP code–level daily preterm birth counts for New York City, 2017–2019 |

## Source

New York City Department of Health and Mental Hygiene. NYSDOH Vital Statistics Birth Micro SAS Datasets. Accessed: 2024-12-13. Available from: https://www.health.ny.gov/statistics/vital_statistics/2019/

## Variables

The file contains ZIP code–level aggregated birth records including:
- ZIP code identifier
- Date of birth
- Preterm birth indicator (gestational age < 37 weeks)
- Total births per ZIP code per day

## Usage

This file is read by `CreateDataFinal.R` (project root) in Step 1 to generate the simulation datasets:
- `truedta_null_Final.csv` — null association scenario (true RR = 1.00)
- `truedta_ass_Final.csv` — harmful association scenario (true RR = 1.115)

## Note

The PM2.5 exposure data (`pm25_js_zip_2000_2016_ny_daily.csv`) is not included in this repository due to its size (117 MB). It is publicly available via Harvard Dataverse: https://doi.org/10.7910/DVN/58C6HG
