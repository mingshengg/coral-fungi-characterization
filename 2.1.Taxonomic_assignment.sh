qiime tools import \
--input-path miseq_mangrove_feb2025.fasta \
--type 'FeatureData[Sequence]' \
--output-path miseq_mangrove_feb2025.qza

#qiime feature-classifier fit-classifier-naive-bayes \
#--i-reference-reads ref_seq.qza \
#--i-reference-taxonomy ref-taxonomy.qza \
#--o-classifier classifier-UNITE-all-eukaryotes-v10.qza

qiime feature-classifier classify-sklearn \
--i-classifier database+classifiers/classifier-UNITE-all-eukaryotes-v10.qza \
--i-reads miseq_mangrove_feb2025.qza \
--o-classification miseq_mangrove_feb2025_taxonomy.qza

qiime tools export \
--input-path rmiseq_mangrove_feb2025_taxonomy.qza \
--output-path exported-taxonomy-shiwen
