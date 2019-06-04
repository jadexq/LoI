//
// Created by WANGJade on 10/12/2018.
//

// steps:
// 0. Generate generate(Constants constants)
// 1. generate.generateX(constants)
// 2. generate.generateY(constants) latent variables
// 3. generate.generateV(constants) manifest variables
// 4. generate.generateXm(constants) cov for missing
// 5. generate.generateR(constants) ind for v_miss

#include "Generate.h"
using namespace std;
using namespace Eigen;

Generate::Generate(Constants constants)
        : x(constants.cn, constants.cp),
          y(constants.cn),
          v(constants.cn, constants.cq),
          r(constants.cn, constants.cq),
          xm1(constants.cn, constants.cr),
          xm(constants.cn, constants.cr+constants.cq+1)
{
    x = MatrixXd::Constant(constants.cn, constants.cp, 0);
    y = VectorXd::Constant(constants.cn, 0);
    v = MatrixXd::Constant(constants.cn, constants.cq, 0);
    r = MatrixXd::Constant(constants.cn, constants.cq, 0);
    xm1 = MatrixXd::Constant(constants.cn, constants.cr, 0);
    xm = MatrixXd::Constant(constants.cn, constants.cr+constants.cq+1, 0);
}

// def generate x
int Generate::generateX(Constants constants)
{
    // simulation
    if (constants.setting == 0)
    {
        VectorXd rnorm_vec = rnorm(constants.cn*constants.cp);
        Map<MatrixXd> rnorm_mat(rnorm_vec.data(), constants.cn, constants.cp);
        // rescale
        ArrayXd rs =  1/(VectorXd::LinSpaced(constants.cp, 1, 5).array());
        x = (rnorm_mat.array().rowwise() * rs.transpose());
    }

    // real data ADNI N=806
    if (constants.setting > 0)
    {
        // read d from file 806*7.txt
//        Eigen::Matrix<double, 806, 7, Eigen::DontAlign> dt;
//        std::ifstream file_dt("d_t.txt"); // d_t.txt 806*7 .txt
//        if (file_dt.is_open())
//        {
//            for (int i = 0; i < 806 * 7; i++)
//                file_dt >> dt(i);
//        }
//        file_dt.close();
        MatrixXd d = MatrixXd::Constant(806,2,0);
        // two tmp
        string str; // store a line as str
        vector<double> vec; // store a line as vec
        // open file
        ifstream file_dt;
        file_dt.open("d.txt");
        // assign
        int i = 0;
        while (getline(file_dt, str))
        {
            // convert str to vec
            strToVec(str, vec);
            // assign row_i in file to col_i in matrix xi
            for (int j=0; j<2; j++)
            {
                d(i,j) = vec[j];
            } // end j

            i++;
        } // end while i
        file_dt.close();

        // read xi from file 806*806, xi_m1.txt, col_i ~ subject_i, xi here row_i ~ subject_i
        MatrixXd xi = MatrixXd::Constant(806, 806, 0);
        // two tmp
        string str2; // store a line as str
        vector<double> vec2; // store a line as vec
        // open file
        ifstream file_xi;
        file_xi.open("xi_z51_mask_70.txt");
        // assign
        i = 0;
        while (getline(file_xi, str2))
        {
            // convert str to vec
            strToVec(str2, vec2);
            // assign row_i in file to col_i in matrix xi
            for (int j=0; j<806; j++)
            {
                xi(j,i) = vec2[j];
            } // end j

            i++;
        } // end while i
        // close file
        file_xi.close();

        // x <- cov | xi, currently no other covariates
        x.rightCols(806) = xi;
        x.leftCols(2) = d;
        // save for checking
        ofstream file;
        file.open("./out/check_xi_read.txt");
        file << xi << endl;
        file.close();

    } // end read d and xi from file
    return 0;
}

// def generate y
// y_i = x_i beta + e_i
int Generate::generateY(Constants constants)
{
    VectorXd rnorm_vec = rnorm(constants.cn);
    VectorXd residual = sqrt(constants.psd) * rnorm_vec;
    VectorXd response = x * constants.beta + residual;
    y = response;

    return 0;
}

// def generate v
// v_i = (1, y_i) clambda + e_i
//     = 1*Lam.row(0) + y_i*Lam.row(1) + e_i
int Generate::generateV(Constants constants)
{
    // simulation
    if (constants.setting == 0)
    {
        VectorXd rnorm_vec = rnorm(constants.cn * constants.cq);
        Map<MatrixXd> rnorm_mat(rnorm_vec.data(), constants.cn, constants.cq);
        v = (rnorm_mat.array().rowwise()*sqrt(constants.psx.transpose().array()) + (y*constants.clambda.row(1)).array()).matrix();
        VectorXd ones = VectorXd::Constant(constants.cn, 1);
        v = (v.array() + (ones * constants.clambda.row(0)).array()).matrix(); // intercept
    }

    // real data ADNI N=806
    if (constants.setting > 0)
    {
        // read v from file vc0_t(treat as continuous) 806*5, keep its shape
        // two tmp
        string str; // store a line as str
        vector<double> vec; // store a line as vec
        // open file
        ifstream file_v;
        file_v.open("v.txt");
        // assign
        int i = 0;
        while (getline(file_v, str))
        {
            // convert str to vec
            strToVec(str, vec);
            // assign row_i in file to col_i in matrix xi
            for (int j=0; j<constants.cq; j++)
            {
                v(i,j) = vec[j];
            } // end j

            i++;
        } // end while i
        file_v.close();

        // save for checking
        ofstream file;
        file.open("./out/check_v_read.txt");
        file << v << endl;
        file.close();
    }

    return 0;
}

// def generateXm
// xm1 n*r, the "ones" are not included
// xm n * (1         + r    + q)
//         intercept   xm1    v
// NOTE, same cov is used for all cq
int Generate::generateXm(Constants constants)
{
    // simulation
    if (constants.setting == 0)
    {
        // xm1, iid N(0,1)
        VectorXd rnorm_vec = rnorm(constants.cn*constants.cr);
        Map<MatrixXd> rnorm_mat(rnorm_vec.data(), constants.cn, constants.cr);
        xm1 = rnorm_mat;
        // col with ones
        VectorXd ones = VectorXd::Constant(constants.cn, 1);
        // combine (ones | xm1 | v) = xm
        xm << ones, xm1, v;
    }

    // real data ADNI N=806
    if (constants.setting > 0)
    {
        // read d from file, 806*7, keep its shape
        MatrixXd d = MatrixXd::Constant(806,2,0);
        // two tmp
        string str; // store a line as str
        vector<double> vec; // store a line as vec
        // open file
        ifstream file_dt;
        file_dt.open("xm.txt");
        // assign
        int i = 0;
        while (getline(file_dt, str))
        {
            // convert str to vec
            strToVec(str, vec);
            // assign row_i in file to col_i in matrix xi
            for (int j=0; j<2; j++)
            {
                d(i,j) = vec[j];
            } // end j

            i++;
        } // end while i
        file_dt.close();

        // save for checking
        ofstream file;
        file.open("./out/check_d_read.txt");
        file << d << endl;
        file.close();

        // col with ones
        VectorXd ones = VectorXd::Constant(constants.cn, 1);

        // combine (ones | xm1 | v) = xm; here xm1 ~ d
        xm << ones, d, v;
    }

    return 0;
}

// def generateR
int Generate::generateR(Constants constants)
{
    // simulation
    if (constants.setting == 0)
    {
        // missing completely at random
        if (constants.model_type==1)
        {
            // binary missing indicator, =1 with prob p_miss
            VectorXd runif_vec = runifVec(constants.cn*constants.cq);
            VectorXd binary_vec = (runif_vec.array()<constants.p_miss).cast<double>();
            Map<MatrixXd> r(binary_vec.data(), constants.cn,constants.cq);
        }

        // missing at random or missing not at randam
        // logit p = xm * betam
        // NOTE under MNAR ind_betam for v all 0
        if (constants.model_type==2|constants.model_type==3)
        {
            // binary missing indicator
            // =1 with prob (1, xm1, v) * (betam0, betam1, betam2)' = xm * betam
            // logit p_miss n*q
            MatrixXd logit_pm = xm * constants.betam;
            // prob miss n*q
            MatrixXd pm = ( exp(logit_pm.array())/(1+exp(logit_pm.array())) );
            // uniform (0,1) matrix n*q
            VectorXd runif_vec = runifVec(constants.cn*constants.cq);
            Map<MatrixXd> runif_mat(runif_vec.data(), constants.cn, constants.cq);
            // r
            r = (runif_mat.array() < pm.array()).cast<double>();
        }
    } // end if simulation

    // real data ADNI N=806
    if (constants.setting > 0)
    {
        // read r from file, 806*5, keep its shape
        // two tmp
        string str; // store a line as str
        vector<double> vec; // store a line as vec
        // open file
        ifstream file_m;
        file_m.open("r.txt");
        // assign
        int i = 0;
        while (getline(file_m, str))
        {
            // convert str to vec
            strToVec(str, vec);
            // assign row_i in file to col_i in matrix xi
            for (int j=0; j<constants.cq; j++)
            {
                r(i,j) = vec[j];
            } // end j

            i++;
        } // end while i
        file_m.close();

        // save for checking
        ofstream file;
        file.open("./out/check_missid_read.txt");
        file << r << endl;
        file.close();
    }

    return 0;
}

// def get x
MatrixXd Generate::getX(void)
{
    return x;
}

// def get y
VectorXd Generate::getY(void)
{
    return y;
}

// def get v
MatrixXd Generate::getV(void)
{
    return v;
}

// def get r, missing indicator of v
MatrixXd Generate::getR(void)
{
    return r;
}

// def get xm, missing cov
MatrixXd Generate::getXm(void)
{
    return xm;
}


// def write x, write x to file x_gen.txt
int Generate::writeX(void)
{
    ofstream file;
    file.open("./out/x_gen.txt");
    file << x << endl;
    file.close();

    return 0;
}

// def write y, write y to file y_gen.txt
int Generate::writeY(void)
{
    ofstream file;
    file.open("./out/y_gen.txt");
    file << y << endl;
    file.close();

    return 0;
}

// def write v, write v to file v_gen.txt
int Generate::writeV(void)
{
    ofstream file;
    file.open("./out/v_gen.txt");
    file << v << endl;
    file.close();

    return 0;
}

// def write r, write r to file r_gen.txt
int Generate::writeR(void)
{
    ofstream file;
    file.open("./out/r_gen.txt");
    file << r << endl;
    file.close();

    return 0;
}

// def write xm, write xm to file xm_gen.txt
int Generate::writeXm(void)
{
    ofstream file;
    file.open("./out/xm_gen.txt");
    file << xm << endl;
    file.close();

    return 0;
}

