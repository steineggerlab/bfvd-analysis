# bfvd-analysis
Scripts to build and analyze the Big Fantastic Virus Database (BFVD), a repository of predicted
structures of viral proteins.

## Versions
| Version | Content | Publication | Code |
|---|---|---|---|
| BFVD v3 | UniProt 2025_03 viral proteins (≤2,000 residues), 5,776,417 structures | to be added | this page |
| BFVD v1 / v2 | UniRef30 2023_02 viral cluster representatives, 351,242 structures | [Nucleic Acids Research, doi:10.1093/nar/gkae1119](https://doi.org/10.1093/nar/gkae1119) | [`script_v1/`](script_v1/README.md) |

## BFVD v3
| Folder | Contents |
|---|---|
| `prediction/` | Data production code for base MSAs, templates, Logan50 search, ColabFold prediction |
| `benchmark/` | benchmark of structure prediction condition |
| `analysis/` | data analysis scripts |
| `common/` | site configuration (`env.example.sh`) |

### Running
```bash
conda env create -f environment.yml && conda activate bfvd-analysis
cp common/env.example.sh common/env.sh                                 # then set paths and databases
```

## Other code
- [`3d-beacons-bfvd/`](3d-beacons-bfvd/): metadata export of BFVD for the 3D-Beacons network.

## Links
- [Webserver](https://bfvd.foldseek.com/): Explore BFVD with UniProt accessions and structures.
- [Data availability](https://bfvd.steineggerlab.workers.dev/): Structures of BFVD and relevant metadata.
