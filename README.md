# Annotate BEDPE Loop Anchors with Genomic Features

**Version:** 2.0  

**Author:** Luis Tenorio Hernández

**Date:** July 2025

---

## Table of Contents

- [Overview](#overview)  
- [Features](#features)  
- [Dependencies](#dependencies)  
- [Usage](#usage)  
- [Input Files](#input-files)  
- [Output](#output)  
- [Workflow Details](#workflow-details)  
- [CLI Options](#cli-options)  
- [Example](#example)  
- [License](#license)

---

## Overview

This script takes a BEDPE file of chromatin loops and annotates each loop anchor with genomic features found in a directory of coordinate files. Supported coordinate formats include:
- BED (e.g., peaks, DHS, ChIP-seq peaks)  
- BEDGRAPH (e.g., signal tracks)
- GTF (to extract gene promoters (TSS +/- 1kb by Default))  

Optionally, in addition to the GTF file you can specify a file with the extension:
- COUNT.CSV (e.g., gene counts in different conditions from RNA-seq)

The result is a TSV of loop anchors with counts of overlaps and associated gene names or signal values per anchor, merged back into BEDPE loop pairs.

---

## Features

- **Anchor extraction** from BEDPE  
- **Intersection with**  
  - BED files → raw overlap counts  
  - BEDGRAPH files → mean signal value per anchor 
  - GTF files → Promoter overlap and gene name aggregation 
- **Gene expression mapping** if a `.counts.csv` file is present alongside GTF inputs  
- **Intermediates** are saved under an `intermediate_files/` folder next to your output, you may delete the after the running.  

---

## Dependencies

- Python 3.7+  
- [pandas](https://pandas.pysample_data.org/)  
- [bedtools](https://bedtools.readthedocs.io/) CLI in your `$PATH`  
- Standard library: `argparse`, `subprocess`, `pathlib`, `io`

---


## Usage

```bash
python3 annotate_loops.py   -b <loops.bedpe>   -d <coords_directory>   -o <annotated_loops.tsv>
```

### CLI Options

| Flag        | Description                                         |
|-------------|-----------------------------------------------------|
| `-b, --bedpe` | Path to input BEDPE file of chromatin loops     |
| `-d, --dir`   | Directory containing coordinate files            |
| `-o, --out`   | Path for the output annotated TSV file           |

---

## Input Files

1. **BEDPE file** (`.bedpe`)  
   - Standard 6-column format:  
     ```
     chr1  start1  end1  chr2  start2  end2  [other cols...]
     ```  
   - Each row defines a chromatin loop; both anchors (A and B) will be annotated.

2. **Coordinate files** in the directory specified by `-d`:
   - **BED** (`*.bed`): genomic intervals to count overlaps.
   - **GTF** (`*.gtf`): gene model annotation; extracts Promoters positions.
   - **BEDGRAPH** (`*.bedgraph`): continuous signal tracks.
   - **Counts CSV** (`*.counts.csv`): table of gene expression counts; used to map expression to overlapped genes for GTF inputs. It's important to ensure your file has this exact extension to be recognized by the script.

---

## Output

- **TSV** (`.tsv`) with one row per loop, columns:
  - `chr_A, start_A, end_A, chr_B, start_B, end_B, loop_id`
  - For each input file:  
    - `<name>_count` for BED counts  
    - `<name>` for comma-separated gene names (GTF)  
    - `<condition>_gene_counts` for expression mapping  
    - `<name>` for mean signal (BEDGRAPH)

- **Intermediate files** saved to:
  ```
  <output_parent>/intermediate_files/
  ├─ anchors.bed
  ├─ <coord_name>_count.bed
  ├─ <coord_name>_Prom.bed
  ├─ <coord_name>_average.bedGraph
  └─ ...
  ```

---

## Workflow Details

1. **Create** a folder for intermediate files `intermediate_files/` beside your output file.  May be erased by the user if not needed
2. **Parse** the BEDPE loops and write a 3-col coordinate BED of all anchors keeping track of the same loop anchors.  
3. **For each** file in `<coords_directory>`:
   - **BED** → count overlaps with `bedtools intersect -c`
   - **GTF** → extract Promoters positions, intersect (`-wa -wb`), group by anchor, list genes, and optionally map expression counts
   - **BEDGRAPH** → intersect (`-wa -wb`), compute mean signal per anchor  
   - **Skip** unsupported extensions with a warning.  
4. **Merge** all annotations back into paired A/B anchors and write the final TSV.

---

## Example

```bash
python Loop_annotation.py   -b sample_loops.bedpe   -d sample_annotation_files  -o sample_loops_annotated.tsv
```

After running, you’ll find:

```text
sample_loops_annotated.tsv
intermediate_files/
├─ anchors.bed
├─ peaks_count.bed
├─ genes_promoter.bed
├─ signal.bedgraph
└─ sample_counts.csv
```

---
