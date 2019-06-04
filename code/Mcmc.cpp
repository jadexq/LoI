//
// Created by WANGJade on 10/12/2018.
//

/*
 how to use this class:
 steps:
 0. Mcmc:mcmc(Constants constants)
 1. readData(Constants constants)
 2. setInitial(Constants constants)
 3. updateBeta(Constants constants)/updatePsi(Constants constants)
 */

#include "Mcmc.h"
using namespace std;
using namespace Eigen;

Mcmc::Mcmc(Constants constants)
        : x(constants.cn, constants.cp),
          y(constants.cn),
          beta(constants.cp),
          ac_beta(constants.cp),
          gamma(constants.cp),
          tau(constants.cp),
          ac_tau(constants.cp),
          a4(constants.cp),
          v(constants.cn, constants.cq),
          r(constants.cn, constants.cq),
          ac_v(constants.cn, constants.cq),
          xm(constants.cn, constants.cr+constants.cq+1),
          betam(constants.cr+constants.cq+1, constants.cq),
          ac_betam(constants.cr+constants.cq+1, constants.cq),
          clambda(2, constants.cq),
          psx(constants.cq)
{
    x = MatrixXd::Constant(constants.cn, constants.cp, 0);
    y = VectorXd::Constant(constants.cn, 0);
    beta = VectorXd::Constant(constants.cp, 0);
    ac_beta = VectorXd::Constant(constants.cp, 0);
    gamma = VectorXd::Constant(constants.cp, 0);
    tau = VectorXd::Constant(constants.cp, 0);
    ac_tau = VectorXd::Constant(constants.cp, 0);
    a4 = VectorXd::Constant(constants.cp, 0);
    v = MatrixXd::Constant(constants.cn, constants.cq, 0);
    r = MatrixXd::Constant(constants.cn, constants.cq, 0);
    ac_v = MatrixXd::Constant(constants.cn, constants.cq, 0);
    xm = MatrixXd::Constant(constants.cn, constants.cr+constants.cq+1, 0);
    betam = MatrixXd::Constant(constants.cr+constants.cq+1, constants.cq, 0);
    ac_betam = MatrixXd::Constant(constants.cr+constants.cq+1, constants.cq, 0);
    clambda = MatrixXd::Constant(2, constants.cq, 0);
    psx = VectorXd::Constant(constants.cq, 0);
}

// def readData, read x/y into Mcmc
int Mcmc::readData(Constants constants, Generate generate)
{
    x = generate.getX();
    //y = generate.getY(); latent variable true value
    v = generate.getV();
    // read missing indicator if there is missing in v
    if (constants.model_type>0)
    {
        r = generate.getR();
        // read cov xm if MAR or MNAR
        if (constants.model_type>1)
        {
            xm = generate.getXm();
        }
    }

    return 0;
}

// def setInitial, set initial value of the para (beta, psd) to be updated
int Mcmc::setInitial(Constants constants, Initial initial)
{
    // prior tau IG(a4, b4=1)
    a4 = initial.sdx;
    a4 = (1/a4.array() + 1).matrix();
    ArrayXd ind = (a4.array() < 3).cast<double>();
    a4 = (3.0*ind + a4.array()*(1-ind)).matrix();
    beta = initial.beta;
    psd = initial.psd;
    gamma = initial.gamma;
    tau = initial.tau;
    theta = initial.theta;
    clambda = initial.clambda;
    psx = initial.psx;
    y = initial.y;
    betam = initial.betam;
    // if there is missing in v, set initial value
    if (constants.model_type>0)
    {
        v = initial.v;
    }
    // if MNAR
    if (constants.model_type==3)
    {
        xm.rightCols(constants.cq) = v;
    }

    return 0;
}

// def getBeta
VectorXd Mcmc::getBeta(void)
{
    return beta;
}

// def getPsd
double Mcmc::getPsd(void)
{
    return psd;
}

// def getGamma
VectorXd Mcmc::getGamma(void)
{
    return gamma;
}

// def getTau
VectorXd Mcmc::getTau(void)
{
    return tau;
}

// def getTheta
double Mcmc::getTheta(void)
{
    return theta;
}

// def getClambda, NOTE clambda must be transposed first
VectorXd Mcmc::getClambda(void)
{
    Map<VectorXd> clam_vec(clambda.data(), clambda.size()); // by col
    return clam_vec;
}

// def getPsx
VectorXd Mcmc::getPsx(void)
{
    return psx;
}

// def getY
VectorXd Mcmc::getY(void)
{
    return y;
}

// def getV
MatrixXd Mcmc::getV(void)
{
    return v;
}

// def getBetam, NOTE betam must be transposed first
VectorXd Mcmc::getBetam(void)
{
    Map<VectorXd> bm_vec(betam.data(), betam.size()); // by col
    return bm_vec;
}


// def updateBeta
// update beta element by element
int Mcmc::updateBeta(Constants constants)
{
    for (int j=0; j<constants.cp; j++)
//    for (int j=0; j<1; j++) // debug
    {
        // copy current beta
        VectorXd beta_p = beta;
        // perturb the jth dim
        // prior 1, random walk, not good
        //beta_p(j) = beta(j) + constants.tun_beta(j)*rnorm(1)(0); // rnorm returns a VectorXd
        // prior 2, half random walk, when gamma_j=0, perturb "0" as the proposal, good
        // beta_p(j) = beta(j)*gamma(j) + constants.tun_beta(j)*rnorm(1)(0); // rnorm returns a VectorXd
        beta_p(j) = beta(j)*gamma(j) + constants.tun_beta(j)*rnorm(1)(0)*sqrt(psd/(x.col(j).array().pow(2).sum()));
        // prior 3, independent MH, use normal to approximate laplace
        // double var = constants.tun_beta(j)*1/( x.col(j).array().pow(2).sum()*(1/psd) + 0.5*gamma(j)*pow(constants.lambda1,2) + 0.5*(1-gamma(j))*pow(constants.lambda0,2) );
        // double mean = 0.25 * var * ((y-x*beta+x.col(j)*beta(j)).array()*x.col(j).array()).sum();
        // cout << "gamj" << gamma(j) << endl;
        // cout << "var" << var << endl;
        // cout << "mean" << mean << endl;
        // beta_p(j) = mean + sqrt(var)*rnorm(1)(0);
        // for calculate loglikelihood value p(Y|all)*p(beta)
        double logli = logPostBeta(beta, constants); // beta
        double logli_p = logPostBeta(beta_p, constants); // beta proposal
        // the acceptance prob
        double prob_ac = exp(logli_p - logli);
        // a U(0,1) random number
        double r01 = runif();
        // ac_ind =1 if accept
        int ac_ind = r01 < prob_ac;
        // update
        if(ac_ind == 1)
        {
            beta(j) = beta_p(j); // if accepted, update beta(j)
            ac_beta(j) = ac_beta(j) + 1; // count the #ac
        }
    }

    return 0;
}

// def updatePsd
int Mcmc::updatePsd(Constants constants)
{
    // the alpha and beta of IG posterior dist
    double ag = constants.a1 + 0.5 * constants.cn;
    double bg = constants.b1 + 0.5 * (y - x * beta).transpose()*(y - x * beta);
    // generate gamma
    double psd_inv = rgamma(1, ag, 1/bg)(0); // length=1 eigen vec
    // inverse
    psd = 1/psd_inv;

    return 0;
}

// def updateGamma
// update gamma element by element
int Mcmc::updateGamma(Constants constants)
{
    for (int j=0; j<constants.cp; j++)
    {
        // calculate the prob of gamma_j = 1
        double nu = constants.lambda1 * exp(-constants.lambda1*abs(beta(j))) * theta;
        double de = nu + constants.lambda0 *tau(j)* exp(-constants.lambda0*tau(j)*abs(beta(j)))*(1-theta);
        double p1 = nu/de;
        // unif(0,1) random number
        double r01 = runif();
        // gam_j = 1 with prob p1
        int gam_j = r01 < p1;
        // assign
        gamma(j) = gam_j;
    }

    // change beta to beta*gamma
    beta = (beta.array() * gamma.array()).matrix();

    return 0;
}

// def updateTau
int Mcmc::updateTau(Constants constants)
{
    // tau proposal
    VectorXd tau_p = tau;
    //tau_p = tau.array() + constants.tun_tau.array()*rnorm(constants.cp).array();
    for (int j=0; j<constants.cp; j++)
    {
        tau_p(j) = tau(j) +  constants.tun_tau(j)*rnorm(1)(0);
        while(tau_p(j) < 0)
        {
            tau_p(j) = tau(j) +  constants.tun_tau(j)*rnorm(1)(0);
        }
    }

    // log likelihood with currently value
    ArrayXd logli = -a4.array()*(tau.array().log()) - constants.lambda0*tau.array()*(beta.array().abs()) - constants.b4/tau.array();
    ArrayXd logli_p = -a4.array()*(tau_p.array().log()) - constants.lambda0*tau_p.array()*(beta.array().abs()) - constants.b4/tau_p.array();
    // ac prob
    ArrayXd prob_ac = exp((logli_p - logli).array());
    // U(0,1) random numbers
    ArrayXd r01 = runifVec(constants.cp);
    // ac_ind =1 if accept
    ArrayXd ac_ind = (r01 < prob_ac).cast<double>();
    // update if gam = 0 & ac_ind =1
    ArrayXd ind = (1-gamma.array()) * ac_ind;
    tau = (tau_p.array()*ind + (1-ind)*tau.array()).matrix();
    // count #ac
    ac_tau = (ac_tau.array() + ind).matrix();

//cout << "taup5 = " << tau_p(6) << endl;
//cout << "tau5 = " << tau(6) << endl;
//cout << "beta5 = " << beta(6) << endl;
//cout << "gam5 = " << gamma(6) << endl;
//cout << "ind5 = " << ind(6) << endl;
//cout << "ac5 = " << ac_tau(6) << endl;

    return 0;
}

// def update Theta, the prob of mixture prior of beta
// posterior Beta(sum_j gam_j + a2, sum_j (1-gam_j) + b2)
int Mcmc::updateTheta(Constants constants)
{
    // calculate ap
    double ap = gamma.sum() + constants.a2;
    // calculate bp
    double bp = (1-gamma.array()).sum() + constants.b2;
    // generate theta|.
    theta = rbeta(1, ap, bp)(0); // rbeta returns a VectorXd

    return 0;
}

// def updateY, the latent variable
// N[mean, sig]
// sig^-1 = 1/psd + clambda * diag(psx_inv) * clambda'
// mean = sig * clambda * diag(psx_inv) * v_i
int Mcmc::updateY(Constants constants)
{
    MatrixXd vcen = (v.array().rowwise() - 1 * clambda.row(0).array()).matrix(); // minus intercept
    double sig = 1 / ( 1/psd + (clambda.row(1).transpose().array().pow(2)*(1/psx.array())).sum() );
    VectorXd mean = sig * ( vcen.array().rowwise() * (clambda.row(1).array()*(1/psx.transpose().array())) ).rowwise().sum();
    mean = mean + (sig/psd) * x * beta;
    VectorXd rnorm_vec = rnorm(constants.cn);
    y = mean + sqrt(sig)*rnorm_vec;

    return 0;
}

// def updateClambda, and psx, the factor loading and residual var of CFA
// v_i = y_i Lambda + e_i
// use joint prior for [Lambda, psx]
// copy WYF's R code "MCMC.R"
int Mcmc::updateClambdaPsx(Constants constants)
{
    // create for later computation
    MatrixXd uk = MatrixXd::Constant(constants.cn, 2, 1); // cov matrix corresponding to free Lam
    VectorXd vkstar = VectorXd::Constant(constants.cn, 0); // the kth dim of v minus fixed r.h.s
    double ap = 0; // 1/psx ~ Gam(ap, bp)
    double bp = 0;
    VectorXd pvar_vec = VectorXd::Constant(2, constants.sig0_clambda); // vec length=2 prior varance sig0_clambda
    MatrixXd pcov_inv = MatrixXd::Constant(2,2,0); // prior cov^-1
    pcov_inv.diagonal() = (1/pvar_vec.array()).matrix();
    MatrixXd cavk = MatrixXd::Constant(2,2,0); // A_vk
    VectorXd avk = VectorXd::Constant(2,0); // a_vk
    double psxk_inv = 0; // 1 / psi_k ~ IG(ap, bp)
    double cavk1 = 0; // when Lam_k fixed (only mu is free), cavk dim = 1
    double avk1 = 0;

    // for Q dim
    for (int k=0; k<constants.cq; k++)
    //for (int k=0; k<1; k++) // debug
    {
        // if Lam_k is free, NOTE the intercept is always free
        if (constants.ind_clambda(1,k)==1)
        {
            vkstar = v.col(k); // because no fixed element
            uk.col(1) = y;
            ap = constants.a3 + 0.5 * constants.cn;
            bp = constants.b3 + 0.5 * vkstar.array().pow(2).sum();
            // A_vk = (H_0vk^-1 + U_k·U_k')^-1, H_0vk is chosen to be diagonal
            // ps: pcov_inv = H_0vk^-1
            cavk = (pcov_inv + uk.transpose()*uk).inverse(); // NOTE .inverse() is better for matrix with dim < 4
            // a_vk = A_vk · (H_0vk^-1·Λ_0vk + U_k·V*_k)
            avk = cavk * ( pcov_inv * constants.clambda0.col(k) + uk.transpose() * v.col(k));
            // β_εk = β_0εk + 0.5(Y*k'·Y*k + Λ_0vk'·H_0vk^-1·Λ_0vk - a_vk'·A_vk^-1·a_vk)
            //      = bp (that we already calculated) + 0.5(Λ_0vk'·H_0vk^-1·Λ_0vk - a_vk'·A_vk^-1·a_vk)
            // ps: β_0εk is b3, in the IG prior for psx
            //     β_εk is bp
            bp = bp + 0.5 * ( constants.clambda0.col(k).transpose()*pcov_inv*constants.clambda0.col(k) - avk.transpose()*cavk.inverse()*avk )(0);
            // ψ_εk|· ~ Inv-Gamma(n/2 + α_0εk, β_εk)
            psxk_inv = rgamma(1, ap, 1/bp)(0); // NOTE our rgamma use a diff parameterization
            // ASSIGN psx(k)
            psx(k) = 1/psxk_inv;
            // Λ_vk|· ~ N[a_vk, ψ_εk·A_vk]
            // ASSIGN clambda(col k)
            clambda.col(k) = mvnorm(1, avk, psx(k) * cavk); //
        }
        else // if Lam_k is fixed, that is only the intercept is free
        {
            vkstar = (v.col(k).array() - clambda(1,k)*y.array()).matrix(); // Lam_k is fixed (at "1")
            ap = constants.a3 + 0.5 * constants.cn;
            bp = constants.b3 + 0.5 * vkstar.array().pow(2).sum();
            // A_vk = (H_0vk^-1 + U_k·U_k')^-1, H_0vk is chosen to be diagonal
            // ps: H_0vk^-1 = 1/sig0_Lam
            //     U_k·U_k' = n
            cavk1 = 1/(1/constants.sig0_clambda + constants.cn);
            // a_vk = A_vk · (H_0vk^-1·Λ_0vk + U_k·Y*_k)
            //wrong-- avk1 = cavk1 * ( (1/constants.sig0_clambda)*constants.clambda0(0,k) + v.col(k).sum() );
            //correct--
            avk1 = cavk1 * ( (1/constants.sig0_clambda)*constants.clambda0(0,k) + vkstar.sum() );
            // β_εk = β_0εk + 0.5(Y*k'·Y*k + Λ_0vk'·H_0vk^-1·Λ_0vk - a_vk'·A_vk^-1·a_vk)
            //      = bp (that we already calculated) + 0.5(Λ_0vk'·H_0vk^-1·Λ_0vk - a_vk'·A_vk^-1·a_vk)
            // ps: β_0εk is b3, in the IG prior for psx
            //     β_εk is bp
            bp = bp+0.5*(constants.clambda0(0,k)*(1/constants.sig0_clambda)*constants.clambda0(0,k)-avk1*(1/cavk1)*avk1);
            // ψ_εk|· ~ Inv-Gamma(n/2 + α_0εk, β_εk)
//cout << "ap=" << ap << " bp=" << bp << endl;
            psxk_inv = rgamma(1, ap, 1/bp)(0); // NOTE rgamma parameterization
            // ASSIGN psx(k)
            psx(k) = 1/psxk_inv;
            // Λ_vk|· ~ N[a_vk, ψ_εk·A_vk]
            // ASSIGN clambda(col k)
            clambda(0,k) = avk1 + sqrt(psx(k)*cavk1) * rnorm(1)(0);
        }
    }

//    cout << "Lam" << endl << clambda << endl;
//    cout << "psx" << endl << psx << endl;

    return 0;
}

// def updateV
// updateV if there is missing in v
int Mcmc::updateV(Constants constants)
{
    // missing completely at randomv OR missing at random
    if (constants.model_type==1|constants.model_type==2)
    {
        // generate v new
        VectorXd rnorm_vec = rnorm(constants.cn * constants.cq);
        Map<MatrixXd> rnorm_mat(rnorm_vec.data(), constants.cn, constants.cq);
        MatrixXd v_new = (rnorm_mat.array().rowwise()*sqrt(psx.transpose().array()) + (y*clambda.row(1)).array()).matrix();
        VectorXd ones = VectorXd::Constant(constants.cn, 1);
        v_new = (v_new.array() + (ones * clambda.row(0)).array()).matrix(); // intercept
        // assign if r = 1
        v = (r.array() * v_new.array() + (1-r.array()) * v.array()).matrix();
//cout << "v: " << v.row(0) << endl;
    }

    // missing not at random
    if (constants.model_type==3)
    {
        // declare for loop use
        ArrayXd expm = ArrayXd::Constant(constants.cn, 0);
        ArrayXd expm_inv = ArrayXd::Constant(constants.cn, 0);
        ArrayXd expm_logit = ArrayXd::Constant(constants.cn, 0);
        ArrayXd pvar = ArrayXd::Constant(constants.cn, 0);
        ArrayXd prop = VectorXd::Constant(constants.cn, 0);
        ArrayXd ll = ArrayXd::Constant(constants.cn, 0);
        ArrayXd llp = ArrayXd::Constant(constants.cn, 0);
        MatrixXd xm_prop = MatrixXd::Constant(constants.cn, 1+constants.cr+constants.cq, 0);
        ArrayXd acp = ArrayXd::Constant(constants.cn, 0);
        ArrayXd r01 = ArrayXd::Constant(constants.cn, 0);
        ArrayXd aci = ArrayXd::Constant(constants.cn, 0);
        VectorXd vj = VectorXd::Constant(constants.cn, 0);

        // for j=1,...,q update v.col(j) for all individual
        // update v.col(j) if r.col(j)=1
        for (int j=0; j<constants.cq; j++)
        {
            // compute proposal var, vec length = n
            // expm = exp( xm * betam.col(j) ), ps: betam dim 1+r+q * q
            expm = exp((xm * betam.col(j)).array()); // vec length = n
            // expm_inv = 1 / (1+expm)
            expm_inv = 1 / (1+expm);
            // expm_logit = expm / (1+expm)
            expm_logit = expm * expm_inv;
            //// proposal var = expm_logit * (expm_inv * betam(1+r+j,j)^2 + 1) - psx(j) WRONG!
            //pvar = expm_logit * (expm_inv * betam(1+constants.cr+j),j + 1) - psx(j);
            // correct: proposal var = 1 / betam(1+r+j,j)^2 * [ expm_logit * (1 - expm_logit) + 1/psx(j) ]
            pvar = 1 / (pow(betam(1+constants.cr+j,j),2) * expm_logit * (1-expm_logit) + 1/psx(j));
            // v.col(j) proposal
            // random walk proposal, add a tun_v(j) to pvar for j=1,...,q
            // prop = v.col(j) + N(0, tun_v(j)*pvar), length=n
            prop = v.col(j).array() + rnorm(constants.cn).array()*(constants.tun_v(j)*pvar.sqrt());
            prop = r.col(j).array()*prop + (1-r.col(j).array())*(v.col(j).array());
//cout << "prop " << j << " : " << endl << prop*(r.col(j).array()) << endl;

            // log likelihood p(v_j miss|.) = p(r_j|.) * P(v_j miss|.)
            // log likelihood of old v.col(j)
            ll = (xm*betam.col(j)).array() - log(1+exp((xm*betam.col(j)).array())) - 0.5*(1/psx(j))*(v.col(j).array()-(clambda(0,j)+y.array()*clambda(1,j))).pow(2);
            // log likelihood of v_j proposal (prop)
            // NOTE when v is changed, xm should also be changed
            xm_prop = xm;
            xm_prop.col(1+constants.cr+j) = prop;
            llp = (xm_prop*betam.col(j)).array() - log(1+exp((xm_prop*betam.col(j)).array())) - 0.5*(1/psx(j))*(prop-(clambda(0,j)+y.array()*clambda(1,j))).pow(2);

            // ac
            acp = exp(llp - ll);
            r01 = runifVec(constants.cn);
            // accept(aci=1) with prob acp
            aci = (r01 < acp).cast<double>();
            // update v.col(j) if: 1. aci=1, 2. r.col(j)=1
            vj = v.col(j); // a copy
            v.col(j) = (prop*(r.col(j).array()*aci) + vj.array()*(1-r.col(j).array()*aci)).matrix();
            // update xm
            xm.col(1+constants.cr+j) = v.col(j);
            // count #ac
            ac_v.col(j) = (ac_v.col(j).array() + r.col(j).array()*aci).matrix();
        }
    }

    return 0;
}

// def updateBetam
// update betam if MAR or MNAR
int Mcmc::updateBetam(Constants constants)
{
    // logit p = xm * betam
    if (constants.model_type==2|constants.model_type==3)
    {
        // declare for loop use
        ArrayXd expm = ArrayXd::Constant(constants.cn, 0);
        ArrayXd expm_inv = ArrayXd::Constant(constants.cn, 0);
        ArrayXd expm_logit = ArrayXd::Constant(constants.cn, 0);
        double pvar = 0;
        double prop = 0; // beta(j,q) proposal
        VectorXd betam_q_proposal = VectorXd::Constant(constants.cr+1, 0);
        double ll = 0;
        double llp = 0;
        double acp = 0;
        double r01 = 0;
        // for q in 1,...,Q (dim v_i = Q)
         for (int q=0; q<constants.cq; q++)
        //for (int q=0; q<1; q++) // debug
        {
            // for r in 1,...,R+Q+1 (dim xm = R+Q+1)
            for (int j=0; j<(constants.cr+constants.cq+1); j++)
            {
                // update betam_j_q if it is free
                if (constants.ind_betam(j,q)==1)
                {
                    // compute proposal variance of betam_j_q
                    // exp(xm_i * betam_q)
                    expm = exp( (xm * betam.col(q)).array() );
                    // 1 / 1+exp(xm_i * betam_q)
                    expm_inv = 1/(1+expm);
                    // exp(xm_i * betam_q) / 1+exp(xm_i * betam_q)
                    expm_logit = expm * expm_inv;
                    // propsal var
                    pvar = 1 / (xm.col(j).array().pow(2)*expm_logit * (1-expm_logit)).sum();

                    // betam_j_q proposal
                    // random walk proposal, add a tun_betam to pvar
                    prop = betam(j,q) + rnorm(1)(0) * sqrt(pvar) * constants.tun_betam(j,q);
                    // create betam_q proposal
                    betam_q_proposal = betam.col(q);
                    betam_q_proposal(j) = prop;

                    // log likelihood p(betam_q|.) = p(r_q|betam_q) * p(betam_q)
                    // log likelihood of old betam_q
                    ll = (r.col(q).array() * (xm*betam.col(q)).array() - log(1+exp((xm*betam.col(q)).array()))).sum();
                    ll = ll - 0.5 * (1/constants.sig0_betam) * pow((betam(j,q)-constants.betam0(j,q)),2);
                    // log likelihood of proposed betam_q
                    llp = (r.col(q).array() * (xm*betam_q_proposal).array() - log(1+exp((xm*betam_q_proposal).array()))).sum();
                    llp = llp - 0.5 * (1/constants.sig0_betam) * pow((prop-constants.betam0(j,q)),2);

                    // ac
                    double acp = exp(llp - ll);
                    double r01 = runif();
                    // accept with prob acp
                    if(r01 < acp)
                    {
                        betam(j,q) = prop; // update betam(j,q)
                        ac_betam(j,q) = ac_betam(j,q) + 1; // count #ac
                    }
                } // end if betam(j,q) is free
            } // end loop j
        } // end loop q
    }

    return 0;
}



// ---------------- private functions ----------------
// calculate the log full conditional dist of Y, log(Y|.)
double Mcmc::logPostBeta(VectorXd b, Constants constants)
{
    // y|.
    double logli1 = -0.5*(1/psd)*(y-x*b).array().pow(2).sum();
    // beta|.
    double logli2 = -((gamma.array()*constants.lambda1+(1-gamma.array())*constants.lambda0*tau.array())*beta.array().abs()).sum();
    double logli = logli1 + logli2;

    return logli;
}

