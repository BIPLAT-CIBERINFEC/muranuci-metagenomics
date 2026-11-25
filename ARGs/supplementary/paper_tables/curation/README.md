# Curated contig-context data

These three CSV files contain the manually selected records and interpretations used to build Supplementary Tables S5A and S5B. They contain public manuscript sample labels (PAC and visit), not the private laboratory-to-patient mapping. They are derived study data, not results calculated automatically by BLAST or by the table generator.

## Files

| File | Unit and contents |
|---|---|
| `curated_contigs.csv` | One row per sample–contig–ARG combination: 11 records representing 10 sample-specific contigs |
| `blast_queries.csv` | One row per reviewed sequence query: 13 records, including contigs and flanking regions |
| `contig_interpretations.csv` | One row per sample-specific contig: 10 records with manually reviewed findings and limitations |

## Contig records

`curated_contigs.csv` contains `sample`, `contig_id`, `contig_length_bp`, `ARG`, `ARG_start`, `ARG_end`, `ARG_strand`, `MMseqs2_assignment` and `MMseqs2_rank`. A contig carrying two ARGs has two rows. Coordinates are one-based, inclusive and refer to the original contig orientation. `ARG_strand` is `+` or `-`. MMseqs2 labels are retained as reported, including GTDB placeholder names.

## Query records and alignment reports

`blast_queries.csv` identifies each query with `sample`, `contig_id`, `BLAST_query_region`, `query_start_in_contig` and `query_end_in_contig`. The generator calculates query length as end minus start plus one. Left and right flanks refer to contig coordinates, irrespective of ARG orientation.

`evidence_available` describes the evidence reviewed. `source_reports` lists any available detailed alignment report filenames, separated by semicolons. Files are read from the sibling `alignment_reports/` directory; they are not included in the public repository. Empty report fields mean that no detailed alignment report is associated with that query here; they do not mean that BLAST found no matches.

`source_report_order` lists the corresponding report order numbers, also separated by semicolons. These retain the row ordering of the reviewed alignments in Table S5B. Multiple reports for one query must have matching filename/order entries. The current data associate nine reports with six complete-contig queries.

## Interpretations

`contig_interpretations.csv` contains `sample`, `contig_id`, `BLAST_and_sequence_comparison_findings` and `interpretation_and_limitations`. These texts record manual assessment and are not generated from hit identity or coverage. Reference-host and plasmid/chromosome descriptions do not establish the bacterial host or genomic location in the metagenomic sample.

## Updating the data

Edit the appropriate CSV to revise a record or interpretation. Preserve sample–contig links across all three files and use UTF-8 encoding. The generator checks required columns, duplicate records, coordinate bounds, interpretation coverage and alignment-report associations. It reads the CSV records in order for S5A and uses `source_report_order` for S5B.

No new database search is performed. Running `generate_tables.py` still requires the private/local analysis inputs and saved alignment reports described in the [ARG analysis README](../../../README.md).
