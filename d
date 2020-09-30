[33m9356b7555ed9f85e887d0a320129a5a5622395b0[m[33m ([m[1;36mHEAD -> [m[1;32mpatch_Section_3.3[m[33m, [m[1;31morigin/patch_Section_3.3[m[33m)[m Vetted ODE_functions
[33m24e95df21c77f133907e4e81d1aeef8aa5e14dbf[m Vetting ODE functions
[33m8bde699069292477483a8f75aecf1f3a7ad564dc[m Vetting -utils- and -eigenvalues-
[33md694ab4b8bc6b59e976841e34748d066567edc62[m Fixed Zero Counting algorithm
[33m3c22920f523b57eee84689e77649046f053c188c[m Previous changes
[33mf2997b17fc1274aa60b9a740e943382ed7519dfb[m Made the eigenvectors symplectically normalized
[33mc19ff09364773b100071c8186a75917df0c51eea[m Removed Parallel processing code
[33m776ce61e771830390747efa0a82938a2b77ce10a[m Tightened up a lot of estimates. Provisionally, it seems that the entire computation can complete successfully for well choosen parameters.
[33m88ce22cca0bfcd1a0f01154a887bd2db071bd1d2[m Added a stable manifold to catch the trajectory
[33m6a097850381ba737e742474b0b40d81996708fec[m Implemented computation putting all of the error in the stable directions
[33m11c9cf6358d20439541eb941c865cd6968da8d09[m[33m ([m[1;31morigin/patch_forward_shooting[m[33m, [m[1;32mpatch_forward_shooting[m[33m)[m Fixed ErrorEigenfunction; added validated eigenvectors
[33m7192f26b4ea71fe6fd3f9f303618ef8c51c62389[m Computing C_G bound from section 3.2
[33me15ca2db0cc8c5ce163de10c3d523d03bafaf01b[m Trying to see if this will fix our problems
[33m5d2aef85f9a7ee257ebbe9736ed3bd8b2e9931fb[m[33m ([m[1;31morigin/feature_CleaningUp_1[m[33m, [m[1;32mfeature_CleaningUp_1[m[33m)[m Cleaned up the Heteroclinic.cpp file
[33m09ada54a30ede9d54b1cd9549361b563814a5755[m Removed unused parts of code
[33m5da831d85be9eb7cdf80ed255673d803c526771b[m[33m ([m[1;31morigin/master[m[33m, [m[1;31morigin/HEAD[m[33m)[m removed unused files
[33m84fa2686da79174771913a3005027876103ddd6f[m .
[33m97ffaa7787b039ea8cf95fe1b687228dbc509654[m This version seems to be fairly stable
[33m47538f24ca27b3564c31efb6f5c87a967f4d936f[m Added plotting of multiple parameters
[33m517210679bb06ca90b926b44dd3e0c72730068fa[m implemented parallel processing
[33mb579cecc20db523248515d7e9c73c2ffce752462[m[33m ([m[1;31morigin/feature_Meet_in_the_Middle_Integration[m[33m)[m Finished this version of multiple shooting, and it works better than single shooting!
[33m6401d45a62762756468a7c167704da3494a6df0d[m[33m ([m[1;31morigin/feature_Vary_Time[m[33m, [m[1;32mfeature_Vary_Time[m[33m)[m It works in n=2, but not n=3
[33mc5090a3f99bf6f58d8ecdaf285a56e380819af74[m I started adding time as a varible to find
[33mf7952be8839605d1f226c050d97c6ee289fca577[m[33m ([m[1;31morigin/patch_Multiple_Shooting_Forward[m[33m, [m[1;32mpatch_Multiple_Shooting_Forward[m[33m)[m Finished changing the multiple shooting so that everything integrates forward
[33m679b18ca2d5988d7f7613c687e994a0f448937da[m Separated out the ODE functions from the main file
[33m4c4590ae3be07b6d08edf9bfc57518c99998aff4[m Added new guesses with the right number of points
[33mfeab7d8b47b2cd4b9ede13deb6f78a8c8099a924[m[33m ([m[1;31morigin/feature_Multiple_Shooting[m[33m, [m[1;32mfeature_Multiple_Shooting[m[33m)[m initial conditions which are verified for n=3
[33m0615bbc0e23f3859454e68d69569184d78fb898a[m Trying to debug the multiple shooting. notwithstanding rate conditions, it can verify n=3, ++
[33m65ad11c196b850bef81bd638b88e51d8eca26a7b[m Multiple shooting is working, not much better than before
[33mc29b806cfab90320e4032ff7452c4d465a103221[m added multiple shooting newton step
[33m78a579c3180076fb72448652191bdba6e32c6ff8[m Continued working on the newton method for multiple shooting
[33md4247c1283f9521c7d20fc2da1a9b1642afe03cc[m Added function to integrate multiple shooting points (IN PROGRESS)
[33mf966de59538ae8471aeb696b20230103e2af1e4a[m Added guess for multiple shooting
[33md60be0eb3b8dfd733a6773f4ed9782aa2888431e[m Created guess for multiple shooting
[33m1acd11d26a4f1644317d71a78856c66dfb6fdea5[m Added an energy function
[33m7d50e0e2af6f7b190abcbf43403c7376eb638162[m[33m ([m[1;32mmaster[m[33m)[m I removed some outputs and tried to parallelize, but then commented out the OpenMP things
[33m5db2a4377bb0091c62a9772f570c7a8d6e9d71d5[m Merge pull request #1 from JCJaquette/JCJaquette-patch-1
[33mcd6d3404bcede24454a3422a3cefda0aad90f273[m[33m ([m[1;31morigin/JCJaquette-patch-1[m[33m)[m Fix Unstable Radius
[33m9f52655ca8c01605014ffbba24384012764a42e3[m Initial Upload
[33m1321fdc312e9b557120951c4b6e9534c18cb0b8b[m Initial commit
