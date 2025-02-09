# Hierarchical Vintage Sparse PCA (Hvsp)

The repository contains code for the "Hierarchical Vintage Sparse PCA (Hvsp): Discussion on Vintage Factor Analysis with Varimax Performs Statistical Inference by Rohe and Zeng". 
It is adapted from [HCD repo](https://github.com/tianxili/HCD).

1. `BTSBM_Hvsp.rmd`: the Rmarkdown file that provides more details on the toy example, simulation and the citation network study. 

2. `BTSBM_Hvsp.html`: knitted HTML of `BTSBM_Hvsp.rmd`

3. `output/`: simulation results. `vsp_2088688.csv` corresponds to the simulation results in 2.1 Varying K in `BTSBM_Hvsp.html`; `vsp_2091690.csv` and `vsp_2099876.csv` correspond to the simulation results in 2.2 Vary degree in `BTSBM_Hvsp.html`

4. `PaperCodeOnline/`: main codes to run simulations. Adapted from [HCD repo](https://github.com/tianxili/HCD).

5. `tables_figures/`: tables and figures generated by `BTSBM_Hvsp.rmd` and used in the discussion.

Data from [HCD repo](https://github.com/tianxili/HCD).

1. `Citation3Core.Rda`: The R data file with the adjacency matrix with row names being author names. It is 707 by 707. It is the pruned core with all nodes have at least three connections, extracted from the largest connected component of the original network

2. `ECV-FullAuthorCommunity.csv`: The full list of authors, their degrees and the cluster label from our algorithm

3. `ECV-Truncated20AuthorCommunity.csv`: The top 20 authors (by degree) from each community (if the community have more than 20 authors) and their research interests from internet