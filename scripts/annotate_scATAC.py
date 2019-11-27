import bio96
import pandas as pd

#load the file containing the location of primer combinations in differend plates:
plate_primer_well = pd.read_csv(snakemake.params[0])

#load the annotation of your experiment:
annotation_df = bio96.load(snakemake.input[0])[['primer_plate','well','plate','cell_type']]

#combine plate with well information to find the correct primer_combination
annotation_df['plate_well'] = annotation_df['primer_plate'] + "_" + annotation_df['well']

#merge dfs to extract primer combi info from plate_primer_well
ann_primers = pd.merge(annotation_df, plate_primer_well, on ='plate_well')
ann_primers['sample_name'] = ann_primers['plate'] + ann_primers['primer_combi']

ann_file = ann_primers[['sample_name','cell_type']]
ann_file.to_csv(snakemake.output[0],index=False)
