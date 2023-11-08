# EIT_Challenge_2023

Authors: 
Daniela Calvetti (dxc57@case.edu), Case Western Reserve University, Cleveland, OH, US <br> 
Monica Pragliola (monica.pragliola@unina.it), University of Naples Federico II, Naples, IT<br> 
Erkki Somersalo (ejs49@case.edu), Case Western Reserve University, Cleveland, OH, US

The codes solve the EIT inverse problem by imposing a sparsity prior on the jumps of the target object. More specifically, the jumps are assumed to follow a conditionally Gaussian distribution with unknown variances, for which a generalized gamma prior is adopted. The resulting hypermodel is minimized by means of the hybrid Iterative Alternating Sequential (IAS) algorithm, originally proposed for linear problems [1] and here extended to the non-linear case.

In line with the Kuopio Tomography Challenge 2023 (KTC2023) guidelines, for which the codes have been designed, the output is the segmentation via the Otsu's method of the reconstructed conductivities.

The main code main_IAS.m takes as input the input folder, the output folder and the difficulty levels (DL), that encodes the downsampling factor in the set of injection patterns. In example.m you can find a call of main_IAS.m for the Training data provided by the KTC2023 organizers. The auxiliary functions and files have to be stored in the same folder as main_IAS.m.




Performance of the proposed algorithm on a target object for different DL.

<img src="https://github.com/MonicaPragliola/EIT_Challenge_2023/assets/122533069/0ca6c35f-d275-450a-940a-a5a52502e02b" width="200">
<img src="https://github.com/MonicaPragliola/EIT_Challenge_2023/assets/122533069/8af11afb-7af8-496b-93f0-c25e86171fe9" width="200">
<img src="https://github.com/MonicaPragliola/EIT_Challenge_2023/assets/122533069/1828433e-60d4-48a3-874b-96106a64081a" width="200">
<img src="https://github.com/MonicaPragliola/EIT_Challenge_2023/assets/122533069/b5028368-a7b3-4b7b-9173-60de37deb419" width="200">




[1] Calvetti D., Pragliola M., Somersalo E. (2020) Sparsity promoting hybrid solvers for hierarchical Bayesian inverse problems, SIAM Journal on Scientific Computing 42(6).



