# Phylogenetic resource composition by Baltimore class

View the [PDF](https://github.com/AngieHinrichs/viral_usher_trees/blob/main/misc/baltimore_classification_plot.pdf)
or [interactive HTML](https://angiehinrichs.github.io/viral_usher_trees/misc/baltimore_classification_plot.html)

baltimore_classification_plot.py uses baltimore_classification.tsv and ../tree_metadata.tsv to plot Baltimore class (plus viroids) vs. tip size of all trees.
To run the script and update baltimore_classification_plot.pdf and baltimore_classification_plot.html:
```bash
pip install -r requirements.txt
python baltimore_classification_plot.py
```

NOTE: when a new tree is added, baltimore_classification.tsv will need to be updated, or else the plots will get an extra column labeled Unclassified.
