This code is associated with the following manuscript:
"Carnelian: alignment-free functional binning and abundance estimation of metagenomic reads"
Sumaiya Nazeen, Bonnie Berger
doi: https://doi.org/10.1101/375121

Upon publication, further information can be found at http://carnelian.csail.mit.edu/

0. Requirments
    Vowpal Wabbit 8.1.1
    scikit-learn
    R 3.3.2
    Python 2.7.13
    BioPython 1.70
    EMBOSS 6.6.0

    This code has been tested with GCC 6.3.0 on Ubuntu 17.04, running
    under Bash 4.4.7(1) on a server with Intel Xeon E5-2695 v2 x86_64 2.40 GHz 
    processor and 320 GB RAM.
    
    Using the EC-2010-DB dataset as gold standard, Carnelian can be comfortably
    run on a machine with 16GB RAM using 1 CPU.

1. Directory structure
data/: EC-2010-DB dataset with gold standard EC labels.
scripts/: R scripts for abundance estimation and analysis from read counts in functional bins.
util/
    ext/: external libararies.
    test/: test drawfrag.c and fasta2skm.c
    drawfrag.c: draw fragments from fasta records.
    fasta2skm.c: construct feature (spaced k-mer profile), and convert to VW input format.
    ldpc.py: generate LSH function using LDPC code.
    sequtil.py: splitting and merging utilities for fasta files.
    merge_pairs.py: links paired-end read files using paired-end relationships.
    reduce.py: translate sequences using reduced amino acid alphabets.
    kseq.h: parse FASTA files
tests/
    demo_data/: data files required for unit and advanced tests.
    basictest_carnelian.py: contains the unit tests for Carnelian.
    advancedtest_carnelian.py: contains the end-to-end tests for Carnelian.
    config.py: configures unit tests and advanced tests for Carnelian.
    README.txt: contains the instructions to run the tests.
    
2. Install and test:
   bash SETUP.sh

3. Usage:

Modes:
    (default --optional-arguments such as k-mer length, fragment size,
    hash functions, etc. are set to work best with EC-2010-DB as used in the manuscript.
    If you're going to train on a different dataset, be sure to tune parameters.)

    1)  ./carnelian.py frag [--optional-arguments] test_dir frag_dir [-h]
        Looks for a fasta file in test_dir with matching label file.
        Randomly draws fragments of length and coverage specified in
        optional-arguments. (use "./carnelian.py frag -h" for details)

        Outputs these fragments with corresponding label into frag_dir.

    2) ./carnelian.py train [--optional-arguments] train_dir model_dir [-h]

        Looks for a fasta file in train_dir with matching label file.
        For each batch of training, randomly draw fragments and generate
        feature vectors using Opal LDPC hashes, and trains Vowpal_Wabbit
        One-Against-All classifier against all batches sequentially. To train 
        classifiers in precise mode, use "--precise" option which will make
        the learned model store probabilities.

        Outputs the generated classifier model into model_dir.

    3) ./carnelian.py retrain [--optional-arguments] old_model_dir new_model_dir new_exmpls_dir [-h]

        Looks for a vowpal-wabbit model with patterns and dictionary file
        in the old_model_dir and a fasta file with matching labels in the
        new_exmpls_dir. Starting with the old model, it updates the existing
        training model and merges new labels with old dictionary using the old
        LDPC patterns. Note that a model trained in default mode must be updated
        in default mode. Same is true for precise mode.
        Output model, dictionary, and pattern files will be generated in new_model_dir.
	
    4) ./carnelian.py translate [--optional-arguments] seq_dir out_dir fgsp_loc [-h]
    	
	Using FragGeneScan program located in the fgsp_loc directory, tries to find
	coding sequences in the input reads fasta file in seq_dir, and translated the
	coding sequences to possible ORFs outputting them in a fasta file in the out_dir.

    5) ./carnelian.py predict [--optional-arguments] model_dir test_dir predict_dir [-h]

        Looks for a classifier model in model_dir, and a fasta file in
        test_dir containing reads/fragments. To make predictions with probabilities,
        run in precise mode using "--precise" option and specify probability cutoff
        using "--cutoff <X>" option.
      
        Outputs the predictions in predict_dir as a fasta file with
        corresponding a corresponding label file.

    6) ./carnelian.py eval reference_file predicted_labels [-h]
    
        Evaluation of prediction accuracy in terms of micro and macro averaged
        precision, sensitivity, and F1-score. If run in "precise" mode, it will
        assume predicted_labels file to have two tab-separated columns: <readID, predLabel>

    7) ./carnelian.py abundance in_dir out_dir mapping_file gs_file [-h]

        Generates abundance estimates of functional terms. Looks for predicted labels for
        each sample in its own sub-directory in in_dir and sample mapping information
        and average protein length per label in mapping_file and gs_file respectively.

        Outputs raw counts and effective counts matrices in out_dir.

    8) ./carnelian.py simulate [--optional-arguments] test_dir train_dir out_dir [-h]

        Runs a full pipeline for performance evaluation starting from training on data 
        in train_dir, testing on data in test_dir, and outputting fragments, model and 
        predictions under out_dir in the following directory structure:

        1frag/
            simulated test data (drawn fragments) are saved here.
            (ignored if --do-not-fragment)
        2model/
            classifier will be saved here.
        3predict/
            fragment classifications are saved here.

    9) ./carnelian.py annotate [--optional-arguments] sample_dir model_dir out_dir fgsp_loc [-h]
	
	Annotates the input nucleotide reads starting from gene finding and translation
	on the reads fasta file in the sample_dir using FragGeneScan located in the fgsp_loc
	directory, then classifying the predicted ORFs using the model in model_dir, and 
	outputting the labels in the out_dir.
	
    Steps to be followed in a typical workflow is given in the workflow.txt file.
    
    To replicate our classification performance analysis the code in performance_analysis.R can be used. Before running
    the script the following packages need to be installed: 
        caret, pROC, ROCR, cvAUC, randomForest

Contact
    Sumaiya Nazeen, nazeen@mit.edu

Acknowledgement
    This implementation of Carnelian is adapted from the source code of the following papers:
    Yunan Luo, Y. William Yu, Jianyang Zeng, Bonnie Berger, and Jian Peng. Metagenomic binning through low density hashing.  Bioinformatics (2018), bty611, https://doi.org/10.1093/bioinformatics/bty611
    K. Vervier, P. Mahe, M. Tournoud, J.-B. Veyrieras, and J.-P. Vert. Large-scale Machine Learning for Metagenomics Sequence Classification , Technical report HAL-01151453, May, 2015.
    
