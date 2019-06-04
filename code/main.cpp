#include <iostream>
#include <ctime>
#include <Eigen/Dense>

#include "Constants.h"
#include "Generate.h"
#include "Initial.h"
#include "Mcmc.h"

using namespace std;
using namespace Eigen;

int main()
{
    // clock
    clock_t start_t, end_t;
    double total_t;
    // starting time
    start_t = clock();


    // -- instantiation "Constants"
    // set: data size, dimensions, hyper-parameters
    Constants constants;
    constants.setIndBetam();

    // store eat of each replication
    // CFA model
    MatrixXd rep_clambda = MatrixXd::Zero(constants.n_rep, 2*constants.cq);
    MatrixXd rep_psx = MatrixXd::Zero(constants.n_rep, constants.cq);
    // LoI
    MatrixXd rep_beta = MatrixXd::Zero(constants.n_rep, constants.cp);
    VectorXd rep_psd = VectorXd::Zero(constants.n_rep);
    MatrixXd rep_gamma = MatrixXd::Zero(constants.n_rep, constants.cp);
    VectorXd rep_theta = VectorXd::Zero(constants.n_rep);
    // missing data mode, for MAR(2) or MNAR(3)
    MatrixXd rep_betam = MatrixXd::Zero(constants.n_rep, constants.cq*(constants.cr+constants.cq+1));

    // replication
    for (int rep=0; rep<constants.n_rep; rep++)
    {
        cout << "-----replication: " << rep+1 << endl;
        // store the posterior samples of para
        // samples_beta(n_mcmc-n_burn, cp)
        MatrixXd samples_beta =  MatrixXd::Constant(constants.n_mcmc-constants.n_burn, constants.cp, 0);
        // samples_psd(n_mcmc-n_burn)
        VectorXd samples_psd =  VectorXd::Constant(constants.n_mcmc-constants.n_burn, 0);
        // samples_gamma(n_mcmc-n_burn)
        MatrixXd samples_gamma =  MatrixXd::Constant(constants.n_mcmc-constants.n_burn, constants.cp, 0);
        // samples_tau(n_mcmc-n_burn)
        MatrixXd samples_tau =  MatrixXd::Constant(constants.n_mcmc-constants.n_burn, constants.cp, 0);
        // samples_theta(n_mcmc-n_burn)
        VectorXd samples_theta =  VectorXd::Constant(constants.n_mcmc-constants.n_burn, 0);
        // samples_clambda(n_mcmc-n_burn, cq*2)
        MatrixXd samples_clambda =  MatrixXd::Constant(constants.n_mcmc-constants.n_burn, 2*constants.cq, 0);
        // samples_psx(n_mcmc-n_burn, cq)
        MatrixXd samples_psx =  MatrixXd::Constant(constants.n_mcmc-constants.n_burn, constants.cq, 0);
        // MAR or MNAR, constants.model_type>1
        cout << "model type --- " << constants.model_type << endl;
        // samples_betam(n_mcmc-n_burn, (cr+cq+1)*cq)
        MatrixXd samples_betam =  MatrixXd::Constant(constants.n_mcmc-constants.n_burn, (constants.cr+constants.cq+1)*constants.cq, 0);

        // -- instantiation "Generate"
        Generate generate(constants);
        // generate data
        generate.generateX(constants);
        generate.generateY(constants);
        generate.generateV(constants);
        // if there is missing in v (manifest variable)
        // MCAR
        if (constants.model_type==1)
        {
            generate.generateR(constants);
        }
        // MAR or MNAR
        if (constants.model_type==2|constants.model_type==3)
        {
            generate.generateXm(constants);
            generate.generateR(constants);
        }

        // write *.txt in ./out
        generate.writeX();
        generate.writeY();
        generate.writeV();
        if (constants.model_type>0)
        {
            generate.writeR();
            if (constants.model_type>1)
            {
                generate.writeXm();
            }
        }

        // -- instantiation "Initial"
        // set initial values of para, betam is always set
        Initial initial(constants);
        // set initial values of y (latent variable)
        initial.setY(constants, generate);
        initial.writeY();
        // if there is missing
        if (constants.model_type>0)
        {
            initial.setV(constants, generate);
            initial.writeV();
        }
        // tau and sdx
        initial.setTauSdx(constants, generate);

        // -- instantiation "Mcmc"
        Mcmc mcmc(constants);
        // read data, given model_type
        mcmc.readData(constants, generate);
        // set initial value of variables to be updated, given model_type
        mcmc.setInitial(constants, initial);

        // update
        for (int g=0; g<constants.n_mcmc; g++)
        {
            // update beta
            mcmc.updateBeta(constants);
            // update psd
            mcmc.updatePsd(constants);
            // update gamma
            mcmc.updateGamma(constants);
            // update tau
            if (constants.ssl == 1)
            {
                mcmc.updateTau(constants);
            }
            // update theta
            mcmc.updateTheta(constants);
            // update y, latent variable
            mcmc.updateY(constants);
            // update clambda and psx in CFA
            mcmc.updateClambdaPsx(constants);
            // if there is missing
            // MCAR
            if (constants.model_type==1)
            {
                // update v
                mcmc.updateV(constants);
            }
            // MAR
            if (constants.model_type==2|constants.model_type==3)
            {
                // update v
                mcmc.updateV(constants);
                // update betam
                mcmc.updateBetam(constants);
            }

            // store the posteror samples
            if (g > (constants.n_burn-1))
            {
                samples_beta.row(g-constants.n_burn) = mcmc.getBeta().transpose();
                samples_psd(g-constants.n_burn) = mcmc.getPsd();
                samples_gamma.row(g-constants.n_burn) = mcmc.getGamma().transpose();
                samples_tau.row(g-constants.n_burn) = mcmc.getTau().transpose();
                samples_theta(g-constants.n_burn) = mcmc.getTheta();
                samples_clambda.row(g-constants.n_burn) = mcmc.getClambda().transpose();
                samples_psx.row(g-constants.n_burn) = mcmc.getPsx().transpose();
                if (constants.model_type>1)
                {
                    samples_betam.row(g - constants.n_burn) = mcmc.getBetam().transpose();
                }
            }
            // print the iterator
            if ((g+1) % 1000 == 0)
            {
                std::cout << "-----iteration: " << g+1 << std::endl;
            }
        } // end loop g

        // print to screen
        VectorXd beta_est = (samples_beta.colwise().sum().array()/(samples_gamma.colwise().sum().array()+1)).matrix();
        double psd_est = samples_psd.mean();
        VectorXd gamma_est = samples_gamma.colwise().mean();
        VectorXd tau_est = samples_tau.colwise().mean();
        double theta_est = samples_theta.mean();
        VectorXd clambda_est_vec = samples_clambda.colwise().mean();
        Map<MatrixXd> clambda_est(clambda_est_vec.data(), 2, constants.cq);
        VectorXd psx_est = samples_psx.colwise().mean();
        std::cout << "beta_est:"<< std::endl << beta_est << std::endl;
        cout << "ac_beta" << std::endl << mcmc.ac_beta.array()/constants.n_mcmc << endl;
        std::cout << "psd_est:" << std::endl << psd_est << std::endl;
        std::cout << "gamma_est:"<< std::endl << gamma_est << std::endl;
        std::cout << "theta_est:" << std::endl << theta_est << std::endl;
        std::cout << "clambda_est:"<< std::endl << clambda_est << std::endl;
        std::cout << "psx_est:"<< std::endl << psx_est << std::endl;
        if (constants.model_type>1)
        {
            VectorXd betam_est_vec = samples_betam.colwise().mean();
            Map<MatrixXd> betam_est(betam_est_vec.data(), constants.cr+constants.cq+1, constants.cq);
            std::cout << "betam_est:"<< std::endl << betam_est << std::endl;
            cout << "ac_betam" << std::endl << mcmc.ac_betam.array()/constants.n_mcmc << endl;
        }
        if (constants.model_type==3)
        {
            cout << "ac_v" << endl << mcmc.ac_v.colwise().sum().array()/(constants.cn*constants.n_mcmc) << endl;
        }
        // Note! can not use burn in when calculate ac_tau
        cout << "ac_tau" << endl << mcmc.ac_tau.array()/(constants.n_mcmc-samples_gamma.colwise().sum().transpose().array());

        // write para mcmc_iter into files
        ofstream file_beta;
        ofstream file_psd;
        ofstream file_gamma;
        ofstream file_tau;
        ofstream file_theta;
        ofstream file_clambda;
        ofstream file_psx;
        file_beta.open("./out/samples_beta.txt");
        file_psd.open("./out/samples_psd.txt");
        file_gamma.open("./out/samples_gamma.txt");
        file_tau.open("./out/samples_tau.txt");
        file_theta.open("./out/samples_theta.txt");
        file_clambda.open("./out/samples_clambda.txt");
        file_psx.open("./out/samples_psx.txt");
        file_beta << samples_beta << endl;
        file_psd << samples_psd << endl;
        file_gamma << samples_gamma << endl;
        file_tau << samples_tau << endl;
        file_theta << samples_theta << endl;
        file_clambda << samples_clambda << endl;
        file_psx << samples_psx << endl;
        file_beta.close();
        file_psd.close();
        file_gamma.close();
        file_tau.close();
        file_theta.close();
        file_clambda.close();
        file_psx.close();
        // if MAR or MNAR
        if (constants.model_type>1)
        {
            ofstream file_betam;
            file_betam.open("./out/samples_betam.txt");
            file_betam << samples_betam << endl;
            file_betam.close();
        }
        // write y (latent variable) est in the last iteration into file
        VectorXd yup = mcmc.getY();
        ofstream file_yup;
        file_yup.open("./out/y_updated.txt");
        file_yup << yup << endl;
        file_yup.close();
        // if there is missing in v, write v updated in the last iteration into file
        if (constants.model_type>0)
        {
            MatrixXd vup = mcmc.getV();
            ofstream file_vup;
            file_vup.open("./out/v_updated.txt");
            file_vup << vup << endl;
            file_vup.close();
        }

        // save the est in this replication
        // CFA
        rep_clambda.row(rep) = clambda_est_vec;
        rep_psx.row(rep) = psx_est;
        // LoI
        rep_beta.row(rep) = beta_est;
        rep_psd(rep) = psd_est;
        rep_gamma.row(rep) = gamma_est;
        rep_theta(rep) = theta_est;
        // missing data model, for MAR(2) and MNAR(3)
        if (constants.model_type > 1)
        {
            VectorXd betam_est_vec2 = samples_betam.colwise().mean(); // "2" avoid redefinition
            rep_betam.row(rep) = betam_est_vec2;
        }
        // write into file
        // CFA
        ofstream file_rep_clambda;
        file_rep_clambda.open("./out/rep_clambda.txt");
        file_rep_clambda << rep_clambda << endl;
        file_rep_clambda.close(); // rep clambda
        ofstream file_rep_psx;
        file_rep_psx.open("./out/rep_psx.txt");
        file_rep_psx << rep_psx << endl;
        file_rep_psx.close(); // rep psx
        // LoI
        ofstream file_rep_beta;
        file_rep_beta.open("./out/rep_beta.txt");
        file_rep_beta << rep_beta << endl;
        file_rep_beta.close(); // rep beta
        ofstream file_rep_psd;
        file_rep_psd.open("./out/rep_psd.txt");
        file_rep_psd << rep_psd << endl;
        file_rep_psd.close(); // rep psd
        ofstream file_rep_gamma;
        file_rep_gamma.open("./out/rep_gamma.txt");
        file_rep_gamma << rep_gamma << endl;
        file_rep_gamma.close(); // rep gamma
        ofstream file_rep_theta;
        file_rep_theta.open("./out/rep_theta.txt");
        file_rep_theta << rep_theta << endl;
        file_rep_theta.close(); // rep theta
        // missing
        ofstream file_rep_betam;
        file_rep_betam.open("./out/rep_betam.txt");
        file_rep_betam << rep_betam << endl;
        file_rep_betam.close(); // rep betam

    } // end loop rep

    // ending time
    end_t = clock();
    total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    std::cout << "total time:" << total_t << std::endl;
    std::cout << "* end *" << std::endl;

    return 0;
}
