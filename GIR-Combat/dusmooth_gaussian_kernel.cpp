#ifndef UTILS_H
#define UTILS_H

#include "Rcpp.h"



void check_subset_vector(Rcpp::IntegerVector, int);

#endif



#include <stdexcept>
#include <sstream>

void check_subset_vector(Rcpp::IntegerVector subvec, int len) {
    for (const auto& s : subvec) {
        if (s==NA_INTEGER || s < 0 || s >= len) {
            throw std::runtime_error("subset indices out of range");
        }
    }
    return;
}


#include "Rcpp.h"



#include <vector>
#include <set>
#include <stdexcept>
#include <cmath>

// [[Rcpp::export]]
SEXP dusmooth_gaussian_kernel(Rcpp::NumericMatrix averaged, Rcpp::IntegerVector index, Rcpp::NumericMatrix mat, double sigma2) {
    const size_t ncells=mat.ncol();
    const size_t ngenes_for_dist=mat.nrow();

    const size_t ngenes=averaged.nrow();
    const size_t nmnn=averaged.ncol();

    if (nmnn!=static_cast<size_t>(index.size())) {
        throw std::runtime_error("'index' must have length equal to number of rows in 'averaged'");
    }


    Rcpp::NumericMatrix output(ngenes, ncells); 
    Rcpp::NumericMatrix exponent(ngenes, ncells);
    std::fill(exponent.begin(), exponent.end(), R_NegInf);

    bool starting_prob=true;
    std::vector<double> distances2(ncells), totalprob(ncells, R_NaReal); 


    auto iIt=index.begin();
    for (size_t i=0; i<nmnn; ++i, ++iIt) {
        auto curcol = mat.column(*iIt);
        auto mnn_iIt = curcol.begin();

        for (size_t other=0; other<ncells; ++other) {
            double& curdist2=(distances2[other]=0);
            double& curdist=(distances2[other]=0);//
            double& curdist2_middle=(distances2[other]=0);//
            auto othercol = mat.column(other);
            auto other_iIt = othercol.begin();
            auto iIt_copy = mnn_iIt;

            for (size_t g=0; g<ngenes_for_dist; ++g) {
                const double tmp=(*iIt_copy  - *other_iIt);
                const double tmp_superman=abs(*iIt_copy - *other_iIt); 
                curdist = std::max(curdist, tmp_superman); 
                curdist2_middle+=tmp*tmp;
                curdist2=0.9*curdist2_middle+0.1*curdist;
                ++other_iIt;
                ++iIt_copy;
            }

           
            curdist2/=-sigma2;
        }
        
      
        double density=0;
        bool starting=true;
        for (const auto& other_mnn : index) {
            if (starting) {
                density=distances2[other_mnn];
                starting=false;
            } else {
                density=R::logspace_add(density, distances2[other_mnn]);
            }
        }

   
        const auto& correction = averaged.column(i);
        auto oIt=output.begin();
        auto eIt=exponent.begin();

        for (size_t other=0; other<ncells; ++other) {
            const double logmult=distances2[other] - density;
            double& logtotal=totalprob[other];

            if (!starting_prob) {
                logtotal=R::logspace_add(logtotal, logmult);
            } else {
                logtotal=logmult;
            }

          
            for (const auto& corval : correction) {
                if (logmult > *eIt) {
                    *oIt *= std::exp(*eIt - logmult);
                    *oIt += corval;
                    *eIt = logmult;
                } else {
                    *oIt += corval * std::exp(logmult - *eIt);
                }
                ++oIt;
                ++eIt;
            }
        }

        starting_prob=false;
    }

 
    for (size_t other=0; other<ncells; ++other) {
        auto curcol=output.column(other);
        auto curexp=exponent.column(other);
        const double logtotal=totalprob[other];

        auto eIt=curexp.begin();
        for (auto& val : curcol) {
            val *= std::exp(*eIt - logtotal);
            ++eIt;
        }
    }

    return output;
}





// #include <Rcpp.h>
// #ifdef RCPP_USE_GLOBAL_ROSTREAM
// Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
// Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
// #endif


// SEXP dusmooth_gaussian_kernel(Rcpp::NumericMatrix averaged, Rcpp::IntegerVector index, Rcpp::NumericMatrix mat, double sigma2);
// RcppExport SEXP sourceCpp_1_dusmooth_gaussian_kernel(SEXP averagedSEXP, SEXP indexSEXP, SEXP matSEXP, SEXP sigma2SEXP) {
// BEGIN_RCPP
//     Rcpp::RObject rcpp_result_gen;
//     Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type averaged(averagedSEXP);
//     Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type index(indexSEXP);
//     Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mat(matSEXP);
//     Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
//     rcpp_result_gen = Rcpp::wrap(dusmooth_gaussian_kernel(averaged, index, mat, sigma2));
//     return rcpp_result_gen;
// END_RCPP
// }
