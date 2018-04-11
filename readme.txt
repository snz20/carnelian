This code is associated with the following manuscript:
"Carnelian: alignment-free functional binning and abundance estimation of metagenomic reads"
Sumaiya Nazeen and Bonnie Berger

Upon publication, further information can be found at http://carnelian.csail.mit.edu/

0. Requirments
    Vowpal Wabbit >= 8.3.1
    scikit-learn
    R 3.3.2

    This code has been tested with GCC 6.3.0 on Ubuntu 17.04, running
    under Bash 4.4.7(1). There are reports of compilation errors on Mac OS X.

1. Directory structure
data/: EC-2192-DB dataset with gold standard EC labels.
util/
    ext/: external libararies.
    test/: test drawfrag.c and fasta2skm.c
    drawfrag.c: draw fragments from fasta records.
    fasta2skm.c: construct feature (spaced k-mer profile), and convert to VW input format.
    ldpc.py: generate LSH function using LDPC code.
    sequtil.py: split a large fasta file into smaller ones and merges after processing.
    kseq.h: parse FASTA files
    2. Install and test:
    bash SETUP.sh

3. Usage: carnelian.py assumes it lives in the current directory structure, but can be symlinked elsewhere.

Modes:
    (default --optional-arguments such as k-mer length, fragment size,
    hash functions, etc. are set to work best with EC-2192-DB as used in the manuscript.
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
        One-Against-All classifier against all batches sequentially. Default
        mode trains SVM classifiers where as precise mode trains logistic
        regression models.

        Outputs the generated classifier model into model_dir.

    3) ./carnelian.py retrain [--optional-arguments] old_model_dir new_model_dir new_examples_dir [-h]

        Looks for a vowpal-wabbit model with patterns and dictionary file
        in the old_model_dir and a fasta file with matching labels in the
        new_examples_dir. Starting with the old model, it updates the existing
        training model and merges new labels with old dictionary using the old
        LDPC patterns. Note that a model trained in default mode must be updated
        in default mode. Same is true for precise mode.
        Output model, dictionary, and pattern files will be generated in new_model_dir.

    4) ./carnelian.py predict [--optional-arguments] model_dir test_dir predict_dir [-h]

        Looks for a classifier model in model_dir, and a fasta file in
        test_dir containing reads/fragments.

        Outputs the predictions in predict_dir as a fasta file with
        corresponding a corresponding label file.

    5) ./carnelian.py eval reference_file predicted_labels [-h]
    
        Evaluation of prediction accuracy in terms of micro and macro averaged
        precision, sensitivity, and F1-score.

    6) ./carnelian.py abundance in_dir out_dir mapping_file gs_file [-h]

        Generates abundance estimates of functional terms. Looks for predicted labels for
        each sample in its own sub-directory in in_dir and sample mapping information
        and average protein length per label in mapping_file and gs_file respectively.

        Outputs raw counts and effective counts matrices in out_dir.

    7) ./carnelian.py simulate [--optional-arguments] test_dir train_dir out_dir [-h]

        Runs a full pipeline training on data in train_dir, testing on
        data in test_dir, and outputting everything under out_dir in the
        following directory structure:

        1frag/
            simulated test data (drawn fragments) are saved here.
            (ignored if --do-not-fragment)
        2model/
            classifier will be saved here.
        3predict/
            fragment classifications are saved here.

Contact
    Sumaiya Nazeen, nazeen@mit.edu

Acknowledgement
    This implementation of Carnelian is adapted from the source code of the following papers:
    Yunan Luo, Y. William Yu, Jianyang Zeng, Bonnie Berger, and Jian Peng. Metagenomic binning through low density hashing. doi: https://doi.org/10.1101/133116
    K. Vervier, P. Mahe, M. Tournoud, J.-B. Veyrieras, and J.-P. Vert. Large-scale Machine Learning for Metagenomics Sequence Classification , Technical report HAL-01151453, May, 2015.
    
