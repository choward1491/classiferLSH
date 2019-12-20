# classiferLSH
This repository represents software and the write up for a final project in CS 598 - Geometric Data Structures under Timothy Chan @ UIUC. This project was basically setup to explore the following:

- Parallel implementations of classical Locality Sensitive Hashing (LSH) and CoveringLSH data structures for approximate r-near neighbor problems using the Hamming metric and l1 metric
- Generalize the r-near neighbor data structures for approximate k-nearest neighbor (k-NN) data structures
- Use the approximate k-NN data structures to build non-parametric classifiers that can return label predictions for some input, given a dataset with points and corresponding labels

The CoveringLSH and classical LSH techniques are covering in [this paper by Pagh in 2016](https://arxiv.org/pdf/1507.03225.pdf). The report for this project is located in the `writeup` directory. The experimental results were run on a 2.3 GHz 8-core Intel Core i9 CPU along with 16 GB 2667 MHz DDR4 RAM and compiled using the latest Clang compiler with support for C++14 and OpenMP.
