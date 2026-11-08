# Metadata snapshot

These tables preserve the pre-consolidation release metadata and updated physical file paths. In the current imaging deposit, group information and complete group-to-animal links are consolidated into `metadata/metadata.csv`; file QC, UID normalization and final radiograph provenance are consolidated into `metadata/files.csv`. Original columns are retained losslessly with group_, qc_, uid_ and provenance_ prefixes; group-animal links use a JSON array field. These snapshot tables support validation provenance, not an alternative animal-identity assignment.
