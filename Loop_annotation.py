# Loop annotation version 2.0 
# This script annotates BEDPE Loop Anchors with Genomic Features.

import argparse  # To parse command-line arguments 
import subprocess  # To execute bash commands (in this script, bedtools)
from pathlib import Path  # To perform directory path operations
import pandas as pd  # To handle dataframes 
import io  # To perform in-memory stream operations


# This function creates a directory to store intermediate files that the user may find useful
# Arguments:
#   base_out: the path to the final output file
# Returns:
#   Path object for the created intermediate_files directory
def make_interm_dir(base_out):
    # Build the path to a new directory named 'intermediate_files'
    interm_dir = Path(base_out).parent / "intermediate_files"
    # Create the directory if it does not already exist
    interm_dir.mkdir(exist_ok=True)
    print(f"Intermediate files will be saved in this directory: {interm_dir}")
    return interm_dir  # Return the path for use by the rest of the workflow


# This function parses a BEDPE file and splits each loop into two anchors,
# labeling them A or B and giving each loop a unique ID.
# Arguments:
#   bedpe_path: path to the BEDPE file
# Returns:
#   DataFrame with columns ['chr','start','end','anchor','loop_id']
def parse_bedpe(bedpe_path):
    print(f"Parsing BEDPE file: {bedpe_path}")
    anchors = []  # Initialize list to hold each anchor tuple
    with open(bedpe_path) as f:
        # Read each line, enumerate gives 1-based loop index
        for i, line in enumerate(f, 1):
            # Split line and take only the first 6 fields (chrom/start/end for both anchors)
            c = line.strip().split()[:6]
            # Append anchor A: chr, start, end, label A, and loop identifier
            anchors.append((c[0], int(c[1]), int(c[2]), "A", f"Loop_{i}"))
            # Append anchor B similarly from columns 4-6
            anchors.append((c[3], int(c[4]), int(c[5]), "B", f"Loop_{i}"))
    # Convert list of tuples into DataFrame for easier downstream handling
    df = pd.DataFrame(anchors, columns=['chr','start','end','anchor','loop_id'])
    return df


# This function writes a DataFrame of anchors to a BED file if it doesn't already exist
# Arguments:
#   df: DataFrame with at least ['chr','start','end'] columns
#   path: Path object for output BED
def write_anchor_bed(df, path):
    if path.exists():
        print(f"Anchor BED already exists: {path}")
    else:
        print(f"Writing anchor BED to: {path}")
        # BED format does not support headers or row indices
        df[['chr','start','end']].to_csv(path, sep='\t', header=False, index=False)


# This function writes a DataFrame of promoters to a BED file if it doesn't exist
# Arguments:
#   prom_df: DataFrame with ['chr','start','end','gene_name']
#   path: Path for output promoter BED
def write_prom_bed(prom_df, path):
    if path.exists():
        print(f"prom BED already exists: {path}")
    else:
        print(f"Writing promoters BED to: {path}")
        # Include gene_name as a 4th column for feature labeling
        prom_df[['chr','start','end','gene_name']].to_csv(path, sep='\t', header=False, index=False)


# This function counts overlaps between anchors and another coordinate file,
# saving counts into a BED file if not already present.
# Arguments:
#   anchor_bed: Path to anchor BED
#   coord: Path to feature BED
#   out: Path to write counts BED
def write_count_bed(anchor_bed, coord, out):
    if out.exists():
        print(f"Count BED already exists: {out}")
    else:
        print(f"Writing count BED to: {out}")
        # Run bedtools intersect with -c to count overlaps per anchor
        subprocess.run([
            'bedtools','intersect','-a',str(anchor_bed),'-b',str(coord),'-c'
        ], check=True, stdout=open(out,'w'))

# This function writes a DataFrame of acerage signal per anchor to a BEDGRAPH file if it doesn't exist
# Arguments:
#   average_df: DataFrame with ['chr','start','end','value']
#   path: Path for saving averages of signal from a BEDGRAPH file
def write_average_bedgraph(average_df, path):
    if path.exists():
        print(f"BEDGRAPH file already exists: {path}")
    else:
        print(f"Writing BEDGRAPH average signal to: {path}")
        # Include value as a 4th column
        average_df[['chr','start','end','value']].to_csv(path, sep='\t', header=False, index=False)



# This function extracts promoter regions from a GTF file,
# defining a +/- 1kb window around transcript start sites.
# Arguments:
#   gtf_path: path to GTF annotation file
#   prom_window: flank size around the TSS to define a promoter, Default: 1000 bp
# Returns:
#   DataFrame with promoter coordinates and gene names
def gtf_to_prom_df(gtf_path, prom_window):
    print(f"Extracting promoters from GTF: {gtf_path}")
    proms = []  # List to hold promoter tuples
    with open(gtf_path) as f:
        for l in f:
            if l.startswith('#'):  # Skip header/comment lines
                continue
            cols = l.rstrip('\n').split('\t')
            # Only keep lines with at least 9 columns and type 'transcript'
            if len(cols) < 9 or cols[2] != 'transcript':
                continue
            attr_field = cols[8]  # The attribute column with gene info
            attrs = {}  # Dictionary to parse key-value pairs
            # Split attributes by semicolon
            for attr in attr_field.split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                parts = attr.split(' ', 1)
                if len(parts) != 2:
                    continue
                key, val = parts
                attrs[key] = val.strip('"')  # Remove surrounding quotes
            # Prefer 'gene_name' if available, otherwise use 'gene_id'
            gene = attrs.get('gene_name') or attrs.get('gene_id')
            try:
                # Determine promoter window based on strand:
                # For '+' strand, promoter initial position is 1kb (Default) upstream of transcript 'start' coordinate
                # For '-' strand, promoter initial position is 1kb (Default) upstream of transcript 'end' coordinate
                if cols[6] == '+':
                    init_pos = max(0, int(cols[3]) - prom_window)
                else:
                    init_pos = max(0,int(cols[4]) - prom_window)
                # Define promoter as +/- 1kb (Default) window
                proms.append((cols[0], init_pos, init_pos + (prom_window*2), gene))
            except ValueError:
                # Skip entries with invalid coordinates
                continue
    # Create DataFrame for promoters
    df = pd.DataFrame(proms, columns=['chr','start','end','gene_name'])
    return df


# This function runs bedtools intersect between two files and returns the result as a DataFrame.
# Arguments:
#   a_path, b_path: file paths for -a and -b inputs
#   flags: list of additional bedtools flags (e.g., ['-wa','-wb'])
#   a_cols: list of column names to assign for A file fields
#   b_cols: list of column names for B file fields
# Returns:
#   DataFrame combining specified columns from both inputs
def run_intersect(a_path, b_path, flags, a_cols, b_cols):
    #print(f"Intersect: -a {a_path} -b {b_path} {' '.join(flags)}")
    try:
        cmd = ['bedtools','intersect','-a',a_path,'-b',b_path] + flags
        # Execute the command and capture stdout for parsing
    except subprocess.CalledProcessError as e:
        SystemError.stderr.write(f"ERROR running bedtools on {b_path.name}:\n{e.stderr}\n")
        exit
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    # Read the output into a DataFrame with no header, assign custom column names
    df = pd.read_csv(io.StringIO(result.stdout), sep='\t', header=None,
                     names=a_cols + b_cols, dtype=str)
    return df


# This function merges the per-anchor annotations back into loop pairs
# Arguments:
#   result: DataFrame with anchor-level annotations including 'anchor' and 'loop_id'
# Returns:
#   DataFrame where each row represents one loop with combined data from both anchors
def make_annotated_paired_df(result):
    # Separate anchor A and anchor B rows
    df_A = result[result['anchor'] == 'A']
    df_B = result[result['anchor'] == 'B']
    # Merge on loop_id to pair A and B annotations side by side
    df_loops = df_A.merge(df_B, on='loop_id', how='inner', suffixes=('_A', '_B'))
    return df_loops


# Main workflow to annotate loop anchors using various genomic feature files
# Arguments:
#   bedpe_file: input BEDPE path
#   coords_dir: directory containing feature files (BED, GTF, bedgraph, counts)
#   output_file: final annotated TSV path
#   prom_window: flank size around the TSS to define a promoter, Default: 1000 bp

def main(bedpe_file, coords_dir, output_file, prom_window=1000):
    print("Starting annotation workflow...")
    # Create directory for intermediate files
    interm_dir = make_interm_dir(output_file)

    # Parse loops into individual anchors
    anchors = parse_bedpe(bedpe_file)
    anchor_bed = interm_dir / 'anchors.bed'
    write_anchor_bed(anchors, anchor_bed)
    # Begin results as copy of anchor-level DataFrame
    results = anchors.copy()

    # Collect any existing counts CSVs into a dictionary for later expression mapping
    counts_files = {p.name: p for p in Path(coords_dir).glob('*.counts.csv')}
    print(f"Counts CSV files: {list(counts_files.keys())}")

    # Loop through each file to intersect in the coordinates directory
    for coord in Path(coords_dir).iterdir():
        # Skip directories
        if not coord.is_file():
            continue
        # Skip COUNTS.CSVs cause they are used only if a gtf file is provided
        if coord.suffix == '.csv' and coord.name.endswith('.counts.csv'):
            continue
        
        # We decompose the name and extension of each file to treat them accordingly
        name, ext = coord.stem, coord.suffix.lower()
        print(f"Processing {coord.name}")

        # If it's a BED feature file, count overlaps per anchor
        if ext == '.bed':
            out = interm_dir / f"{name}_count.bed"
            write_count_bed(anchor_bed, coord, out)
            # Read counts and add as new column in results
            df = pd.read_csv(
                out, sep='\t', header=None,
                names=['chr','start','end','count'], usecols=['chr','start','end','count'],
                dtype={'count': int}
            )
            results[name] = df['count']

        # If it's a GTF annotation, extract promoters and intersect
        elif ext == '.gtf':
            prom_df = gtf_to_prom_df(coord, prom_window)
            prom_bed = interm_dir / f"{name}_prom.bed"
            write_prom_bed(prom_df, prom_bed)

            # Intersect anchors with promoter BED, keeping both sets of columns
            df = run_intersect(
                str(anchor_bed), str(prom_bed), ['-wa','-wb'],
                ['chr','start','end'], ['prom_chr','prom_start','prom_end','gene_name']
            )
            # Drop duplicates in case the same gene is hit multiple times by one anchor
            df = df.drop_duplicates(subset=['chr', 'start', 'end', 'gene_name'])

            # Convert types
            df['start'] = df['start'].astype(int)
            df['end'] = df['end'].astype(int)

            # Merge to results — this will create multiple rows per anchor if multiple genes match
            results = results.merge(
                df[['chr', 'start', 'end', 'gene_name']],
                on=['chr', 'start', 'end'],
                how='left'
            )

            # Rename the column to the current feature name
            results.rename(columns={'gene_name': name}, inplace=True)

            # If expression count files exist, map counts to each gene
            if counts_files:
                csv_path = next(iter(counts_files.values()))
                print(f"Mapping counts from {csv_path.name}")
                counts_df = pd.read_csv(csv_path)
                
            # Identify gene column and sample columns
            gene_col = counts_df.columns[0]
            cond_cols = list(counts_df.columns[1:])
            
            # Merge all expression values at once
            results = results.merge(counts_df, left_on=name, right_on=gene_col, how='left')

            # Rename columns to indicate they are gene counts
            for cond in cond_cols:
                results.rename(columns={cond: f"{cond}_gene_counts"}, inplace=True)

            # Optionally drop gene_col (original from counts_df)
            results.drop(columns=[gene_col], inplace=True)

        # If it's a bedgraph file, compute mean signal over each anchor
        elif ext == '.bedgraph':
            # Intersect anchors with bedgraph to collect value per overlap
            df = run_intersect(
                str(anchor_bed), str(coord), ['-wa','-wb'],
                ['chr','start','end'], ['bg_chr','bg_start','bg_end','value']
            )
            # Convert value column to numeric for averaging
            df['value'] = pd.to_numeric(df['value'], errors='coerce')
            grouped = df.groupby(['chr','start','end'])['value'].mean().reset_index()
            grouped['start'] = grouped['start'].astype(int)
            grouped['end'] = grouped['end'].astype(int)
            grouped[name] = grouped['value'].fillna(0)
            bedgraph_file = interm_dir / f"{name}_average.bedGraph"
            write_average_bedgraph(grouped, bedgraph_file)
            # Merge mean signal into results
            results = results.merge(
                grouped[['chr','start','end', name]], on=['chr','start','end'], how='left'
            )

        else:
            # Skip any unsupported file types
            print(f"Skipping unsupported file: {coord.name}")

    # After processing all feature files, reconstruct per loop annotation
    df_loops = make_annotated_paired_df(results)
    print(f"Writing output to {output_file}")
    # Write loop table to TSV, including both anchors' annotations
    df_loops.to_csv(output_file, sep='\t', index=False)
    print("Workflow completed.")


# --- Command-line interface ---
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate bedpe loop anchors with intersection features')
    parser.add_argument('-b', '--bedpe', required=True, help='Path to input BEDPE file')
    parser.add_argument('-d', '--dir', required=True, help='Directory with coordinate files')
    parser.add_argument('-o', '--out', required=True, help='Output annotated file (tsv)')
    parser.add_argument('-pw', '--promwindow', required=False, help='Flank size, in base pairs, around the TSS to define promoters. Default: +/- 1000')
    args = parser.parse_args()
    main(args.bedpe, args.dir, args.out)
