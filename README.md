## Extracting Brain Disease Related Connectome Subgraphs by Adaptive Dense Graph Discovery

This is an implementation of ADSD algorithm to identify connectomic subgraphs that are systematically associated with brain disorders from group-level brain connectome data.

# Codes

"greedy_ad.m" is the function for greedy algorithm (Algorithm (1) in Web Appendix) to maximize the adaptive density function with a given lambda.
"greedy_ad_blik.m" is the complete ADSD algorithm with lambda selected by integrated likelihood (Algorithm (2) in Web Appendix). 
"sim_clique_repeat_wgold.m" is for replicated simulation results with comparison to Charikar's and Goldberg's algorithm.

# Examples
We include a small simulation example in "demo.m". This file includes the generating mechanism of a simulation data set and the implementation of ADSD. To see the demonstrative document with codes and output, please refer to "demo.html" observed from the publish option of Matlab. "W1_correct.mat" includes the edge-wise inference results for a brain network data from schizophrenia study
