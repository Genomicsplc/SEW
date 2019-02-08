// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
using namespace Rcpp;



//' @export
// [[Rcpp::export]]
Rcpp::List cpp_single_sample_expectation(
    const Rcpp::List& sampleReads,
    const arma::mat& eHapsCurrent,
    const int nSNPs,
    const int N_r,
    const bool calculate_updates = true,
    const bool calculate_read_probabilities = true,
    int t_break = 0,
    int t_reading_offset = 0,
    int t_writing_offset = 0
) {
    // do not allow K to vary in this model
    const int K = 2;
    //
    // initialize output. make stub output if not needed depending on updating
    //
    double log_likelihood = 0;
    int N_r_output = 1;    
    if (calculate_read_probabilities) {
        N_r_output = N_r;
    }
    int nSNPs_output = 1;
    if (calculate_updates) {
        nSNPs_output = nSNPs;
    }
    arma::rowvec p_reads = arma::zeros(1, N_r_output);
    arma::mat p_reads_given_hap_k = arma::zeros(N_r_output, 2);
    arma::mat p_h_given_O = arma::zeros(N_r_output, 2);
    arma::mat eHapsUpdate_numer = arma::zeros(nSNPs_output,K);
    arma::mat eHapsUpdate_denom = arma::zeros(nSNPs_output,K);
    //
    // initialize intermediate variables
    //
    int k, iRead, j, J, u_j, use_k, u_j_reading, u_j_writing;
    double pR, pA, p_read, eps, p_r_no_j, p_h_and_g_1, p_h_and_g_2;
    arma::rowvec p_read_given_hap_k = arma::zeros(1, K);    
    //
    // loop
    //
    for(int iRead=0; iRead<= N_r - 1; iRead++) {
        Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
        J = as<int>(readData[0]); // 0-based number of Unique SNPs on read
        arma::ivec bqU = as<arma::ivec>(readData[2]); // bq for each SNP
        arma::ivec u = as<arma::ivec>(readData[3]); // position of each SNP from 0 to nSNPs-1
        arma::mat pr = arma::zeros(J + 1, 2);
        for(j=0; j<=J; j++) {
            if(bqU(j)<0) {
                eps = pow(10,(double(bqU(j))/10));
                pR=1-eps;
                pA=eps/3;
            }
            if(bqU(j)>0) {
                eps = pow(10,(-double(bqU(j))/10));
                pR=eps/3;
                pA=1-eps;
            }
            pr(j,0)=pR; // do not need to re-initialize - J controls that
            pr(j,1)=pA;
        }
        // flip_k is a vector with 0 for no flip and 1 for a flip
        arma::rowvec flip_k = arma::zeros(1, J + 1);
        if (t_break > 0) {
            for(j=0; j<=J; j++) {
                if (t_break <= u(j))
                    flip_k(j)=1;
            }
        }
        arma::mat p_read_given_hap_k_components = arma::zeros(J + 1, K);        
        for(j=0; j <= J; j++) {
            u_j_reading = u(j) - t_reading_offset;            
            for(k=0; k < K; k++) {
                use_k = k;
                if (flip_k(j) == 1)
                    use_k = 1 - k; // 0-based; k=0 -> 1, k=1 -> 0
                p_read_given_hap_k_components(j, k) = \
                    eHapsCurrent(u_j_reading, use_k) * pr(j,1) + \
                    (1 - eHapsCurrent(u_j_reading, use_k)) * pr(j,0);
            }
        }
        p_read_given_hap_k.fill(1);
        for(j=0; j <= J; j++) {
            for(k=0; k < K; k++) {
                p_read_given_hap_k(k)=p_read_given_hap_k(k) * \
                    p_read_given_hap_k_components(j, k);
            }
        }
        
        p_read = 0.5 * (p_read_given_hap_k(0) + p_read_given_hap_k(1));
        log_likelihood = log_likelihood + -log(p_read);

        // save outputs
        if (calculate_read_probabilities) {
            p_reads(iRead)=p_read;
            p_reads_given_hap_k.row(iRead)=p_read_given_hap_k;
        }
        // do updates
        if (calculate_updates) {
            for(j=0; j <= J; j++) {
                u_j = u(j);
                u_j_reading = u(j) - t_reading_offset;
                u_j_writing = u(j) - t_writing_offset;
                for(k=0; k < K; k++) {
                    use_k = k;
                    if (flip_k(j) == 1)
                        use_k = 1 - k; // 0-based; k=0 -> 1, k=1 -> 0
                    p_r_no_j = p_read_given_hap_k(k) / p_read_given_hap_k_components(j,k);
                    p_h_and_g_1 = p_r_no_j * 0.5 * pr(j, 1) *   \
                        eHapsCurrent(u_j_reading, use_k) / p_read;
                    p_h_and_g_2 = p_r_no_j * 0.5 * pr(j, 0) *   \
                        (1 - eHapsCurrent(u_j_reading, use_k)) / p_read;
                    if (calculate_read_probabilities) {
                        p_h_given_O(iRead, k) = p_h_and_g_1 + p_h_and_g_2;
                    }
                    eHapsUpdate_numer(u_j_writing, k) = eHapsUpdate_numer(u_j_writing, k) + \
                        p_h_and_g_1;
                    eHapsUpdate_denom(u_j_writing, k) = eHapsUpdate_denom(u_j_writing, k) + \
                        p_h_and_g_1 + p_h_and_g_2;
                }
            }
        }

    }

    //
    // return
    //
    if (calculate_read_probabilities == false & calculate_updates == true) {
        return(wrap(Rcpp::List::create(
            Rcpp::Named("eHapsUpdate_numer") = eHapsUpdate_numer,
            Rcpp::Named("eHapsUpdate_denom") = eHapsUpdate_denom,
            Rcpp::Named("log_likelihood") = log_likelihood
        )));
    } else if (calculate_read_probabilities == true & calculate_updates == true) {
        return(wrap(Rcpp::List::create(
            Rcpp::Named("eHapsUpdate_numer") = eHapsUpdate_numer,
            Rcpp::Named("eHapsUpdate_denom") = eHapsUpdate_denom,
            Rcpp::Named("log_likelihood") = log_likelihood,            
            Rcpp::Named("p_reads") = p_reads,
            Rcpp::Named("p_reads_given_hap_k") = p_reads_given_hap_k,
            Rcpp::Named("p_h_given_O") = p_h_given_O
        )));
    } else if (calculate_updates == false) {
        return(wrap(Rcpp::List::create(
            Rcpp::Named("log_likelihood") = log_likelihood
        )));
    }
    // not in use, prevents warnings I think
    return(wrap(Rcpp::List::create(
        Rcpp::Named("log_likelihood") = log_likelihood
    )));
}
