library(qs)
library(SeuratObject)
library(Matrix)

WORKDIR <- '/Users/jeffrey.katz/Documents/side-projects/bioml/mice-exercise'
TMP <- file.path(WORKDIR, 'seurat_tmp')
dir.create(TMP, showWarnings = FALSE)

cat('Reading .qs file...\n')
obj <- qread(file.path(WORKDIR, 'GSE183288_Single_cell_atlas.qs'))
cat('Done reading.\n')

rna <- obj@assays[['RNA']]

# Write counts matrix
cat('Writing counts matrix...\n')
writeMM(rna@counts, file.path(TMP, 'counts.mtx'))
writeLines(rownames(rna@counts), file.path(TMP, 'genes.txt'))
writeLines(colnames(rna@counts), file.path(TMP, 'barcodes.txt'))

# Write normalized data if different from counts
if (!identical(rna@data, rna@counts) && length(rna@data) > 0) {
  cat('Writing normalized data matrix...\n')
  writeMM(rna@data, file.path(TMP, 'data.mtx'))
}

# Write metadata
cat('Writing metadata...\n')
write.csv(obj@meta.data, file.path(TMP, 'metadata.csv'))

# Write reductions
for (red_name in names(obj@reductions)) {
  red <- obj@reductions[[red_name]]
  cat('Writing reduction:', red_name, '\n')
  write.csv(red@cell.embeddings, file.path(TMP, paste0(red_name, '.csv')))
}

cat('All data exported to', TMP, '\n')
