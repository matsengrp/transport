from collections import Counter
import os
from random import sample
from statsmodels.distributions.empirical_distribution import ECDF

import numpy as np
import pandas as pd

from config import CONFIG
from python.tcr_scorer import TCRScorer

class RandomizationTest():
    output_dir = CONFIG["TMP_OUTPUT"]
    trial_1_file = os.path.join(output_dir, "trial_1.csv")
    trial_2_file = os.path.join(output_dir, "trial_2.csv")
    
    def __init__(self, df_1, df_2, observed_scores, species, trial_count=100):
        self.df_1 = df_1
        self.df_2 = df_2
        self.observed_scores = observed_scores
        self.species = species
        self.N1 = df_1.shape[0]
        self.N2 = df_2.shape[0]
        self.df_2_tcrs = list(dict.fromkeys(self.df_2['tcr']))
        self.full_df = pd.concat([self.df_1, self.df_2], axis=0).reset_index(drop=True)
        self.all_indices = self.full_df.index.tolist()
        self.trial_count = trial_count
        self.do_randomization_test()    

    def split_datasets(self):
        df_1_trial_indices = sample(self.all_indices, self.N1)
        df_1_trial = self.full_df.iloc[df_1_trial_indices]
        df_2_trial = self.full_df.iloc[~self.full_df.index.isin(df_1_trial_indices)]
    
        df_1_trial.iloc[:, [0, 1]].to_csv(self.trial_1_file, header=False, index=False)
        df_2_trial.iloc[:, [0, 1]].to_csv(self.trial_2_file, header=False, index=False)
    
        return df_1_trial, df_2_trial
    
    def do_randomization_test(self):
        self.enrichment_dict = {tcr: [] for tcr in self.df_2_tcrs}
    
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        for trial in range(self.trial_count):
    
            df_1_trial, df_2_trial = self.split_datasets()
    
            trial_scorer = TCRScorer(
                file_1=self.trial_1_file,
                file_2=self.trial_2_file,
                species=self.species
            )

            os.remove(self.trial_1_file)
            os.remove(self.trial_2_file)
    
            for tcr, score in trial_scorer.enrichment_dict.items():
                if tcr in self.df_2_tcrs:
                    self.enrichment_dict[tcr].append(score)

        min_sample_size = np.min([len(x) for x in self.enrichment_dict.values()])

        for tcr, scores in self.enrichment_dict.items():
            downsampled_scores = sample(scores, min_sample_size)
            tcr_ecdf = ECDF(downsampled_scores)
            tcr_obs_score = self.observed_scores[tcr]
            self.enrichment_dict[tcr] = {
                "scores": downsampled_scores,
                "z_score": (tcr_obs_score - np.mean(downsampled_scores))/np.std(downsampled_scores),
                "p_value": 1 - tcr_ecdf(tcr_obs_score),
            }
