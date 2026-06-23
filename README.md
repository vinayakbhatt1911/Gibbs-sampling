# Gibbs-sampling
DNA Motif Discovery using Gibbs Sampling
Project Overview
This project implements the Gibbs Sampling algorithm for discovering conserved DNA motifs across multiple DNA sequences.
The algorithm iteratively refines motifs by:
Building a profile matrix with pseudocounts
Computing k-mer probabilities
Sampling motifs probabilistically
Updating motifs to improve motif conservation
Concepts Used
DNA Sequences
k-mers
Profile Matrix (PWM)
Pseudocounts
Gibbs Sampling
Motif Scoring
Probabilistic Search
Functions
random_kmer()
Selects a random k-mer from a DNA sequence.
build_profile()
Constructs a profile matrix using pseudocounts.
get_prob()
Computes the probability of a k-mer using the profile matrix.
profile_random_kmer()
Selects a new motif according to weighted probabilities.
find_score()
Calculates motif conservation score.
Gibbs()
Main Gibbs Sampling algorithm.