<!-- Modeled after the Tidyverse news file -->
<!-- https://style.tidyverse.org/news.html -->

# loomR 0.2.1

## LOOM_SPEC_VERSION

 - 2.0.1

## Breaking changes

 - Removed support for `loom$batch.scan` and `loom$batch.next`, chunk generation is now handled by `loom$chunk.points` and `loom$chunk.indices`
 - Changed parameter names to better match Seurat v3 (see section "New parameter names below")

## New Features

### New parameter names

 - loomR 0.2.1 brings about several new parameter names to better match Seurat v3
 
| Old Name | New Name |
| -------- | -------- |
| `display.progress` | `verbose` |
| `dataset.use` | `dataset` |
| `genes.use` | `features` |
| `cells.use` | `cells` |
| `gene.names` | `feature.names` |
| `gene.attrs` | `feature.attrs` |
| `do.transpose` | `transpose` |
| `layer` | `layers` |
| `attribute` | `attributes` |
| `attribute.names` | `attributes`|
| `calc.numi` | `calc.count` |


## Minor improvements and fixes

 - Removed dependencies on iterators and itertools
 - Change 'nGene' to 'nFeatures' and 'nUMI' to 'nCount' when calculating nUMI/nCount and nGene/nFeatures upon object creation
 - Changed default cell names attribute to 'CellID'
 - Changed default gene names attribute to 'Gene'
