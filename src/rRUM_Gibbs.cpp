#include <RcppArmadillo.h>
#include <rgen.h>
#include <simcdm.h>

// [[Rcpp::interfaces(r, cpp)]]

arma::field<arma::mat>
parm_updatecpp(unsigned int N, unsigned int J, unsigned int K, unsigned int C,
               const arma::mat Y, const arma::mat &Q, arma::mat &alpha,
               arma::cube &X, arma::mat &Smat, arma::mat &Gmat, arma::vec &pi,
               const arma::vec vv, const arma::vec &delta0, double as = 1,
               double bs = 1, double ag = 1, double bg = 1)
{

    arma::vec k_index = arma::linspace(0, K - 1, K);
    double kj, prodXijk, pi_ijk, aik, u, compare;
    double pi_ik, aik_nmrtr_k, aik_dnmntr_k, c_aik_1, c_aik_0;
    arma::vec aik_nmrtr(K);
    arma::vec aik_dnmntr(K);

    // Update X cube
    for (unsigned int i = 0; i < N; ++i) {

        arma::vec ai = (alpha.row(i)).t();
        arma::vec Yi = (Y.row(i)).t();
        arma::vec ui = arma::randu<arma::vec>(K);
        aik_nmrtr = arma::ones<arma::vec>(K);
        aik_dnmntr = arma::ones<arma::vec>(K);

        for (unsigned int j = 0; j < J; ++j) {

            double Yij = Yi(j);
            arma::vec Xij = X.tube(i, j);
            arma::uvec task_ij = find(Q.row(j) == 1);

            for (unsigned int k = 0; k < task_ij.n_elem; ++k) {
                kj = task_ij(k);
                aik = alpha(i, kj);
                Xij(kj) = 1;
                prodXijk = prod(Xij(task_ij));
                u = R::runif(0.0, 1.0);
                pi_ijk = (1.0 - prodXijk) * (aik * (1.0 - Smat(j, kj)) +
                                             (1.0 - aik) * Gmat(j, kj));
                compare = (pi_ijk > u);
                Xij(kj) = (1.0 - Yij) * compare + Yij;

                aik_nmrtr(kj) = (Xij(kj) * (1.0 - Smat(j, kj)) +
                                 (1.0 - Xij(kj)) * Smat(j, kj)) *
                                aik_nmrtr(kj);
                aik_dnmntr(kj) = (Xij(kj) * Gmat(j, kj) +
                                  (1.0 - Xij(kj)) * (1.0 - Gmat(j, kj))) *
                                 aik_dnmntr(kj);
            }
            X.tube(i, j) = Xij;
        }

        // Update alpha_ik
        for (unsigned int k = 0; k < K; ++k) {
            ai(k) = 1.0;
            c_aik_1 = (arma::conv_to<double>::from(ai.t() * vv));
            ai(k) = 0.0;
            c_aik_0 = (arma::conv_to<double>::from(ai.t() * vv));

            aik_nmrtr_k = aik_nmrtr(k) * pi(c_aik_1);
            aik_dnmntr_k = aik_dnmntr(k) * pi(c_aik_0);
            pi_ik = aik_nmrtr_k / (aik_nmrtr_k + aik_dnmntr_k);
            ai(k) = 1.0 * (pi_ik > ui(k));
        }
        alpha.row(i) = ai.t();
    }

    // update pi
    arma::vec a_bijection = alpha * vv;
    arma::uvec deltatilde =
        arma::hist(a_bijection, arma::linspace<arma::vec>(0, C - 1, C));
    pi = rgen::rdirichlet(deltatilde + delta0);

    // update Smat and Gmat
    arma::vec pistar = arma::zeros<arma::vec>(J);
    arma::mat rstar = arma::zeros<arma::mat>(J, K);
    double pg, ps, ug, us, gjk, sjk;

    for (unsigned int j = 0; j < J; ++j) {
        arma::uvec task_ij = find(Q.row(j) == 1);
        arma::mat Xj = X.tube(0, j, N - 1, j);
        double pistar_temp = 1.0;

        for (unsigned int k = 0; k < task_ij.n_elem; ++k) {
            kj = task_ij(k);
            arma::vec Xjk = Xj.col(kj);
            arma::vec ak = alpha.col(kj);

            double Sumalphak = (arma::conv_to<double>::from(ak.t() * ak));
            double SumXjk = (arma::conv_to<double>::from(Xjk.t() * Xjk));
            double SumXjkalphak = (arma::conv_to<double>::from(Xjk.t() * ak));
            double bsk = SumXjkalphak;
            double ask = Sumalphak - SumXjkalphak;
            double agk = SumXjk - SumXjkalphak;
            double bgk = N - SumXjk - Sumalphak + SumXjkalphak;
            ug = R::runif(0.0, 1.0);
            us = R::runif(0.0, 1.0);

            // draw g conditoned upon s_t-1
            pg = R::pbeta(1.0 - Smat(j, kj), agk + ag, bgk + bg, 1, 0);
            gjk = R::qbeta(ug * pg, agk + ag, bgk + bg, 1, 0);
            // draw s conditoned upon g
            ps = R::pbeta(1.0 - gjk, ask + as, bsk + bs, 1, 0);
            sjk = R::qbeta(us * ps, ask + as, bsk + bs, 1, 0);

            Gmat(j, kj) = gjk;
            Smat(j, kj) = sjk;

            rstar(j, kj) = gjk / (1.0 - sjk);        // compute rstarjk
            pistar_temp = (1.0 - sjk) * pistar_temp; // compute pistarj
        }
        pistar(j) = pistar_temp;
    }

    arma::field<arma::mat> out(4);
    out(0) = pistar, out(1) = rstar, out(2) = pi, out(3) = alpha;

    return out;
}

Rcpp::List rrum_main(const arma::mat &Y, const arma::mat &Q,
                     const arma::vec &delta0, unsigned int chain_length = 10000,
                     double as = 1, double bs = 1, double ag = 1, double bg = 1)
{

    unsigned int N = Y.n_rows;
    unsigned int J = Y.n_cols;
    unsigned int K = Q.n_cols;
    unsigned int C = static_cast<unsigned int>(pow(2.0, static_cast<double>(K)));

    arma::vec vv = simcdm::attribute_bijection(K);

    // Prior values for betas and Dirichlet distribution
    // arma::vec delta0 = arma::ones<arma::vec>(C);

    // Savinging output
    arma::mat PISTAR(J, chain_length);
    arma::cube RSTAR(J, K, chain_length);
    arma::mat PIs(C, chain_length);
    arma::cube ALPHAS(N, K, chain_length);

    // need to initialize, alphas, X,ss, gs,pis
    arma::mat alpha = arma::randu<arma::mat>(N, K); // K>1 is assumed
    alpha.elem(find(alpha > 0.5)).ones();
    alpha.elem(find(alpha <= 0.5)).zeros();
    arma::cube X = arma::zeros<arma::cube>(N, J, K);
    arma::mat ss = arma::randu<arma::mat>(J, K);
    arma::mat gs =
        (arma::ones<arma::mat>(J, K) - ss) % arma::randu<arma::mat>(J, K);
    arma::vec pis = rgen::rdirichlet(delta0);

    // Start Markov chain
    for (unsigned int t = 0; t < chain_length; ++t) {
        // updata X,alpha,pi,s,g,pistar,rstar
        arma::field<arma::mat> output =
            parm_updatecpp(N, J, K, C, Y, Q, alpha, X, ss, gs, pis, vv, delta0);

        // update value for pis. alphas are updated via pointer. save classes
        // and PIs
        PISTAR.col(t) = output(0);
        RSTAR.slice(t) = output(1);
        PIs.col(t) = output(2);
        ALPHAS.slice(t) = output(3);
    }

    return Rcpp::List::create(
        Rcpp::Named("PISTAR", PISTAR), Rcpp::Named("RSTAR", RSTAR),
        Rcpp::Named("PI", PIs), Rcpp::Named("ALPHA", ALPHAS));
}

//' Gibbs sampler to estimate the rRUM
//'
//' Obtains samples from posterior distributon for the reduced Reparametrized
//' Unified Model (rRUM).
//'
//' @param Y            A `matrix` with \eqn{N} rows and \eqn{J} columns
//'                     indicating the indviduals' responses to each of the
//'                     items.
//' @param Q            A `matrix` with \eqn{J} rows and \eqn{K} columns
//'                     indicating which attributes are required to answer each
//'                     of the items.An entry of 1 indicates attribute \eqn{k}
//'                     is required to answer item \eqn{j}. An entry of one
//'                     indicates attribute \eqn{k} is not required.
//' @param chain_length A `numeric` indicating the number of iterations of
//'                     Gibbs sampler to be run.  Default is set to 10000.
//' @param as           A `numeric`, parameter for the prior distribution of
//'                     `pistar`.  High values as encourage higher values of
//'                     `pistar` and lower values of `rstar`.
//' @param bs           A `numeric`, parameter for the prior distribution of
//'                     `pistar`.  High values as encourage lower values of
//'                     `pistar` and higher values of `rstar`.
//' @param ag           A `numeric`, parameter for the prior distribution of
//'                     `rstar`.  High values as encourage higher values of
//'                     `rstar`.
//' @param bg           A `numeric`, parameter for the prior distribution of
//'                     `pistar`.  High values as encourage lower values of
//'                     `rstar`.
//' @param deltas       A `vector`, parameters for the Dirichlet prior on `pi`.
//'
//' @return A `List`
//'
//' - `PISTAR` A `matrix` where each column represents one draw from the
//'            posterior distribution of pistar.
//' - `RSTAR` A \eqn{J x K x chain_length} `array` where `J` reperesents the
//'           number of items, and `K` represents the number of attributes.
//'           Each slice represents one draw from the posterior distribution
//'           of `rstar`.
//' - `PI` `matrix` where each column reperesents one draw from the posterior
//'        distribution of `pi`.
//' - `ALPHA` An \eqn{N x K x chain_length} `array` where `N` reperesents the
//'           number of individuals, and `K` represents the number of
//'           attributes. Each slice represents one draw from the posterior
//'           distribution of `alpha`.
//' @author 
//' Steven Andrew Culpepper, Aaron Hudson, and James Joseph Balamuta
//' 
//' @template rrum-example
//' @template rrum-references
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List rrum_helper(const arma::mat &Y, const arma::mat &Q,
                       const arma::vec &delta0,
                       unsigned int chain_length = 10000, double as = 1,
                       double bs = 1, double ag = 1, double bg = 1)
{

    // Consider improving error messages... 
    if (Q.n_rows != Y.n_cols) {
        Rcpp::stop("`Y` must have as many rows as `Q` has columns");
    }
    
    if (delta0.n_elem !=
        static_cast<unsigned int>(pow(2.0, static_cast<double>(Q.n_cols)))) {
        Rcpp::stop("`delta0` must be numeric of length 2 ^ ncol(Q)");
    }

    return(rrum_main(Y, Q, delta0, chain_length, as, bs, ag, bg));
}
