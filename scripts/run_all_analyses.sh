rm -rf output
rm -rf tmp_output

echo "Running replicate analysis..."
python analyses/replicates.py
Rscript --vanilla --slave R/cutoff.R

echo "Running z_score analysis..."
python analyses/z_scores.py
Rscript --vanilla --slave R/z_score.R

echo "Running motif analysis..."
python analyses/motif.py
Rscript --vanilla --slave R/motif.R
