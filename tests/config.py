import os

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"demo_data")

small_fasta_file = os.path.join(data_dir, "small.fasta")
small_fasta_num_records = 100
small_label_file = os.path.join(data_dir, "small.label")
small_label_pred_file = os.path.join(data_dir, "pred.label")
small_label_prob_pred_file = os.path.join(data_dir, "small-prob.label")
small_label_num_labels = 100
small_label_unique_labels = 5

old_dico_file1 = os.path.join(data_dir, "a-dico.txt")
old_dico_file2 = os.path.join(data_dir, "c-dico.txt")
new_dico_file1 = os.path.join(data_dir, "b-dico.txt")
new_dico_file2 = os.path.join(data_dir, "d-dico.txt")
merged_dico_file = os.path.join(data_dir, "merged-dico.txt")

two_col_file = os.path.join(data_dir, "two-col.tab")
one_col_file = os.path.join(data_dir, "one-col.txt")

small_vw_file = os.path.join(data_dir, "small.vw")
small_vw_prob_file = os.path.join(data_dir, "prob.vw")
small_label_prob_file = os.path.join(data_dir, "prob.label")
small_label_prob_tsv_file = os.path.join(data_dir, "label-probabilities.tsv")
cutoff = 0.80

dummy_aa_in_dir = os.path.join(data_dir, "aa_in")
dummy_aa_out_dir = os.path.join(data_dir, "aa_out")
dummy_gs_file = os.path.join(data_dir, "gs_data.tab")
dummy_mapping_file = os.path.join(data_dir, "dummy_map.tab")

orf_seq = os.path.join(data_dir, "orf.fasta")
orf_label = os.path.join(data_dir, "orf.label")
orf_precise_label = os.path.join(data_dir, "orf_precise.label")
read_seq = os.path.join(data_dir, "reads.fasta")
reads_label = os.path.join(data_dir, "reads.label")
reads_precise_label = os.path.join(data_dir, "reads_precise.label")
reads_ref_label = os.path.join(data_dir, "reads_wrt_ref.label")
reads_ref_label_precise = os.path.join(data_dir, "reads_wrt_ref_precise.label")

one_fasta = os.path.join(data_dir, "one.fasta")
one_peptide_fasta = os.path.join(data_dir, "one.peptides.fasta")

splits_dir = os.path.join(data_dir, "splits")
ftest_dir = os.path.join(data_dir, "frag_test")
frag_dir = os.path.join(data_dir, "frag_test/tmp")
ttest_dir = os.path.join(data_dir, "trans_test")
pep_dir = os.path.join(data_dir, "trans_test/tmp")

model_dir = os.path.join(data_dir, "model")
model_precise_dir = os.path.join(data_dir, "model_precise")

predict_dir = os.path.join(data_dir, "predictions")
predict_precise_dir = os.path.join(data_dir, "predictions_precise")

new_examples_dir = os.path.join(data_dir, "new_examples")
new_model_dir = os.path.join(data_dir, "new_model")
new_model_precise_dir = os.path.join(data_dir, "new_model_precise")

samples_dir = os.path.join(data_dir,"samples")
map_file = os.path.join(data_dir,"map.tab")
