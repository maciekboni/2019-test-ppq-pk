# 2019-test-ppq-pk

Pharmacokinetic-Pharmacodynamic model of anti-malarial drug resistance to *Plasmodium falciparum*

## Lumefantrine pk-pd model based on Kloprogge et al., 2018
https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002579

## Usage:
 `run_ppq_pk [drug module] [options]`

Example: `run_ppq_pk --AL -n 100`


### Available Drug Modules

`--adq`     : Amodiaquine

`--art`     : Artemisinin monotherapy

`--AL`      : Artemether-Lumefantrine 

`--lum`     : Lumefantrine

### Options

`-o`        : Output format. Default outputs hourly drug concentration for all patients, as well as            parasitemia.

              1 - Output only the pmax and corresponding efficacies.

`-n`        : Number of patients (default: 1)

`--age`     : Age of patients

`--weight`  : Weight of patients

`--pmf`     : Parasite Multiplication Factor (default: 10.0 for a 48h cycle)


`--hill_art`: Hill co-efficient for Dihydroartemisinin in the Artemether-Lumefantrine module

`--ec50_art`: EC50 for Dihydroartemisin in the Artemether-Lumefantrine module

`--pmax_art`: PMax for Dihydroartemisin in the Artemether-Lumefantrine module


`--hill_lum`: Hill co-efficient for Lumefantrine in the Artemether-Lumefantrine module

`--ec50_lum`: EC50 for Lumefantrine in the Artemether-Lumefantrine module

`--pmax_lum`: PMax for Lumefantrine in the Artemether-Lumefantrine module