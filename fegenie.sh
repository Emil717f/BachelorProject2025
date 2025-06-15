conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict


conda init
conda create -n fegenie -c conda-forge -c bioconda -c defaults fegenie=1.0 --yes
conda activate fegenie
FeGenie.py -h

conda activate fegenie

FeGenie.py -bin_dir test_genomes -bin_ext fna


cp */*.fna 


# Run this
FeGenie.py -bin_dir all_fna_genomes -bin_ext fna # Runs FeGenie on every fna file in the all_fna_genomes folder, -which contains extracted fna files from all_genomes/extracted_files folder



