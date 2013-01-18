I. Quick start.

The merge_metadata.py script has been included in the MaAsLin package to help add metadata to otu tables (or any tab delimited file where columns are the samples). This script was used to make the maaslin_demo.pcl file found in this project.

The generic command to run the merge_metadata.py is:
python merge_metadata.py input_metadata_file < input_measurements_file > output_pcl_file

An example of the expected files are found in this project in the directory maaslin/input/for_merge_metadata
An example of how to run the command on the example files is as follows (when in the maaslin folder in a terminal):
python src/merge_metadata.py input/for_merge_metadata/maaslin_demo_metadata.metadata < input/for_merge_metadata/maaslin_demo_measurements.pcl > input/maaslin_demo.pcl

II. Script overview
merge_metadata.py takes a tab delimited metadata file and adds it to a otu table. Both files have expected formats given below. Additionally, if a pipe-delimited consensus lineage is given in the IDs of the OTUs (for instance for the genus Bifidobacterium, Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Bifidobacterium), the higher level clades in the consensus lineage are added to other otu in the same clade level generating all higher level clade information captured in the otu data*. This heirarchy is then normalized using the same heirarchical structure. This means, after using the script, a sample will sum to more than 1, typically somewhere around 6 but will depend on if your data is originally at genus, species, or another level of resolution. All terminal otus (or the original otus) in a sample should sum to 1. 

*To help combat multiple comparisons, additional clades are only added if they add information to the data set. This means if you have an otu Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Bifidobacterium and no other related otus until Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales, Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae will not be added to the data set because it will be no different than the already existing and more specific Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Bifidobacterium otu. Clades at and above Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales will be included depending on if there are other otus to add them to at those clade levels.


III. Description of input files

Metadata file:
Please make the file as follows:
1. Tab delimited
2. Rows are samples, columns are metadata
3. Sample Ids in the metadata file should match the sample ids in the otu table.
4. Use NA for values which are not recorded.
5. An example file is found at input/for_merge_metadata/maaslin_demo_metadata.metadata

OTU table:
Please make the file as follows:
1. Tab delimited.
2. Rows are otus, columns are samples (note this is transposed in comparison to the metadata file).
3. If a consensus lineage is included in the otu name, use pipes as the delimiter.
4. An example file is found at input/for_merge_metadata/maaslin_demo_measurements.pcl
