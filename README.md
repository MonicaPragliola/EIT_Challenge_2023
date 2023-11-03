# EIT_Challenge_2023

Authors: Daniela Calvetti (dxc57@case.edu), Case Western Reserve University, Cleveland, OH, US 
         Monica Pragliola (monica.pragliola@unina.it), University of Naples Federico II, Naples, IT
         Erkki Somersalo (ejs49@case.edu), Case Western Reserve University, Cleveland, OH, US

The codes solve the EIT inverse problem by imposing a sparsity prior on the jumps of the target object. More specifically, the jumps are assumed to follow a conditionally Gaussian distribution with unknown variances, for which a generalized gamma prior is adopted. The resulting hypermodel is minimized by means of the hybrid Iterative Alternating Sequential (IAS) algorithm, originally proposed for linear problems [1] and here extended to the non-linear case.

In line with the Kuopio Tomography Challenge 2023 (KTC2023) guidelines, for which the codes have been designed, the output is the segmentation via the Otsu's method of the reconstructed conductivities.

The main code main.m runs the function IAS_for_EIT.m on the four training data provided by the KTC2023 committee for different acquisition set-up. The quality of the segmented output conductivities is assessed by means of a scoring function designed by the KTC2023 committee.



[1] Calvetti D., Pragliola M., Somersalo E. (2020) Sparsity promoting hybrid solvers for hierarchical Bayesian inverse problems, SIAM Journal on Scientific Computing 42(6).


Performance of the proposed algorithm on a target object for different difficulty levels (DL), i.e. for different and decreasing numbers of injection patterns.

<img src="https://github.com/MonicaPragliola/EIT_Challenge_2023/assets/122533069/0ca6c35f-d275-450a-940a-a5a52502e02b" width="220">
<img src="https://github.com/MonicaPragliola/EIT_Challenge_2023/assets/122533069/8af11afb-7af8-496b-93f0-c25e86171fe9" width="220">
<img src="https://github.com/MonicaPragliola/EIT_Challenge_2023/assets/122533069/1828433e-60d4-48a3-874b-96106a64081a" width="220">
<img src="https://github.com/MonicaPragliola/EIT_Challenge_2023/assets/122533069/b5028368-a7b3-4b7b-9173-60de37deb419" width="220">






