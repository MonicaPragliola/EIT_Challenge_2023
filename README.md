# EIT_Challenge_2023

Authors: Daniela Calvetti (dxc57@case.edu), Case Western Reserve University, Cleveland, OH, US 
         Monica Pragliola (monica.pragliola@unina.it), University of Naples Federico II, Naples, IT
         Erkki Somersalo (ejs49@case.edu), Case Western Reserve University, Cleveland, OH, US

The codes solve the EIT inverse problem by imposing a sparsity prior on the jumps of the target object. More specifically, the jumps are assumed to follow a conditionally Gaussian distribution with unknown variances, for which a generalized gamma prior is adopted. The resulting hypermodel is minimized by means of the hybrid Iterative Alternating Sequential (IAS) algorithm, originally proposed for linear problems [1] and here extended to the non-linear case.

In line with the Kuopio Tomography Challenge 2023 (KTC2023) guidelines, for which the codes have been designed, the output is the segmentation via the Otsu's method of the reconstructed conductivities.

The repository contains:
- 'TrainingData' folder, containing the data (provided by the KTC2023 committee)
- 'GroundTruths' folder, containing the reference segmentations (provided by the KTC2023 committee)
- 'Output' folder, in which the segmented conductivities are saved
- 'MiscCodes' folder, containing a miscellaneous of auxiliary codes, ranging from the segmentation algorithm to the scoring function evaluating the results (provided by the KTC2023 committee)
- auxiliary files EITFEM.m, sigmaplotter.m, Mesh_dense.mat, Mesh_sparse.mat (provided by the KTC2023 committee)
- auxiliary files GetEdges.m, GetIncrementMatrix.m, CGLS_update_x0
- the Matlab function IAS_for_EIT.m
- the Matlab script main.m, that runs IAS_for_EIT.m on the four target image on different level of difficulties, i.e. for different acquisition set-up

The user can simply run the main script main.m.






[1] Calvetti D., Pragliola M., Somersalo E. (2020) Sparsity promoting hybrid solvers for hierarchical Bayesian inverse problems, SIAM Journal on Scientific Computing 42(6).
