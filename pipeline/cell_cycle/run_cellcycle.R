library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))
source(here('pipeline', 'cell_cycle', 'assign_cell_cycle.R'))

parser <- ArgumentParser(description = "Assign cell cycle phases")

parser$add_argument('--sce_clus', metavar='DIRECTORY', type='character', help="path to sce_clus")
parser$add_argument('--cell_cycle_species', metavar='FILE', type='character', help='human or mouse sample: supply "human" or "mouse"')
parser$add_argument('--output_file_name', metavar='FILE', type='character', help='where we save the sce object with cell cycle gene information')
parser$add_argument('--corrected_file_name', metavar='FILE', type='character', help='where we save the sce object after correcting for cell cycle effect')

args <- parser$parse_args()


assign_cell_cycle(path_to_sce=args$sce_clus, 
                  species=args$cell_cycle_species, 
                  output_path=args$output_file_name,
                  corrected_path=args$corrected_file_name)

