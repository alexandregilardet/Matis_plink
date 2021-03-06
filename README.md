# Matis_plink

## Description

plinkPy takes .ped and .map files exported from Plink as input. After parsing of both files into a single object in an OOP-manner, 
it allows to extract part of the data and export it back into .ped and .map files as well as csv, genepop and structure formats.

## Getting started

When the `plinkPy.py` file is in the working directory with an empty `__init__.py` file:
```
from plinkPy import *
```
### Input
The .ped .map input should display A/B alleles with 0 as the missing value.

## Usage

### Create the plinkPy instance and parse the files
```
plink_1 = plinkPy('plink_export.ped', 'plink_export.map')
```

### Access properties
Get the individuals as well as their number:
```
plink_1.individuals
```
Get the markers as well as their number:
```
plink_1.markers
```

### Filter
By a list of individuals `target_1`:
```
individuals_of_interest = plink_1.get_individuals(target_1)
```
By a list of markers `target_2`:
```
markers_of_interest = plink_1.get_markers(target_2)
```

### Export
```
individuals_of_interest.export(format, 'filtered_export')
```
The `format` can be:
- 'ped_map'
- 'csv'
- 'genepop'
- 'structure'

Export to csv, genepop and structure changes the original order of markers into an alphabetic one.

Export to genepop and structure adds an empty blank line at the end of the file.

Export to genepop misses a blank space after the comma at individual lines.

Export to structure codes missing value as `9`.

