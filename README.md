The Source Data and the code for reproducing the figures in TaxVAMB paper.

`data` folder contains all the necessary Source Data.

Most of the figures of the paper are made in Python and can be reproduced by running the Jupyter notebook `scripts_python/Figures.ipynb`.
```
pip install -r requirements.txt
jupyter notebook
```

The trees from Figure 4c, Figure 6a and Supplementary Figure 8 are made in R. The Figure 4c for instance only requires the summary and the tree that is an output of GTDBtk (v220 here, so `gtdbtk.bac120.summary.tsv` and `gtdbtk.backbone.bac120.classify.tree`).
```
Rscript scripts_R/Figure4c.R --tree data/human_longread.tree --summary data/human_longread_summary.tsv --output figures/figure4c_human_gut.pdf
Rscript scripts_R/Figure4c.R --tree data/sludge.tree --summary data/sludge_summary.tsv --output figures/figure4c_sludge.pdf
Rscript scripts_R/Figure6a.R
Rscript scripts_R/Supplementary_Figure8.R
```

The figures will appear in the `figures` folder.