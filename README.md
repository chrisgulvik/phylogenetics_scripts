# Phylogenetics Scripts

##### System Requirements
- biopython
- PhyML
- RAxML compiled with the PTHREADS lib

##### Example Installation
    pip install biopython
    brew install phyml raxml
    git clone https://github.com/chrisgulvik/phylogenetics_scripts.git $HOME/phylogenetics_scripts
    echo 'export PATH="$PATH:$HOME/phylogenetics_scripts"' >> $HOME/.bash_profile

##### Usage
- Build a ML tree using __make_ML_tree.py__

    `make_ML_tree.py -f input.fasta.aln -e .fasta.aln -T 15`
