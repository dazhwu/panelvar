//
// Created by Tiger on 7/11/2022.
//

#include "IRF.h"


RowMatrixXd PAR1_matrix(RowMatrixXd &beta, int lags) {
  // cols(): num of dep     rows(): num of indep

	int num_dep = beta.cols();
	int num_indep = num_dep * lags;

	////cout << num_indep << " " << num_dep << endl;

  RowMatrixXd tbp = RowMatrixXd::Zero(num_indep, num_indep);

  if (lags == 1) {
    tbp.block(0, 0, num_indep, num_dep) = beta.block(0, 0, num_indep, num_dep);

  } else {
    tbp.block(0, 0, num_indep, num_dep) = beta.block(0, 0, num_indep, num_dep);
    for (int i = 0; i < num_indep - num_dep; ++i) {
      tbp(i, i + num_dep) = 1;
    }
  }

  return (tbp);
}

vector<RowMatrixXd> IRF_matrix(RowMatrixXd &beta, int ahead, int lags){

	int num_dep = beta.cols();
	int num_indep = num_dep * lags;

	RowMatrixXd iden=RowMatrixXd::Identity(num_indep,num_indep);
	RowMatrixXd par1= PAR1_matrix(beta, lags);
	RowMatrixXd J=RowMatrixXd::Zero(num_dep, num_indep);
	for(int i =0;i<num_dep;++i)
		J(i,i)=1;
	vector<RowMatrixXd> tbr;
	tbr.reserve(ahead);
	//MatrixPower<RowMatrixXd> Apow(par1);

	tbr.push_back(RowMatrixXd::Identity(num_dep, num_dep));

	RowMatrixXd temp= iden;
	for(int i=1; i<ahead; ++i){
		temp=temp* par1;
		RowMatrixXd new_mat=J *temp* J.transpose();
		tbr.push_back(new_mat);

	}


	return tbr;

}

RowMatrixXd vcov_residual(RowMatrixXd &residual, VectorXi &is_na ){
	int num_obs=is_na.size()-is_na.sum();
	RowMatrixXd nona_residual(num_obs, residual.cols());
	int new_row=0;
	for(int i=0;i<residual.rows(); ++i){
		if (is_na(i)==0){

			nona_residual.row(new_row)=residual.row(i);
			new_row+=1;
		}
	}

	RowMatrixXd centered = nona_residual.rowwise() - nona_residual.colwise().mean();
	RowMatrixXd cov = (centered.adjoint() * centered) / double(nona_residual.rows() - 1);

	return cov;
}

vector<RowMatrixXd>irf(string method, RowMatrixXd &residual, VectorXi &is_na, RowMatrixXd &beta, int ahead, int lags){
	if (method=="girf")
		return girf(residual, is_na, beta,  ahead, lags);
	else
		return oirf(residual, is_na, beta,  ahead, lags);
}

vector<RowMatrixXd>oirf(RowMatrixXd &residual, VectorXi &is_na, RowMatrixXd &beta, int ahead, int lags){
	vector<RowMatrixXd> ma_phi=IRF_matrix(beta, ahead, lags);

	RowMatrixXd cov= vcov_residual(residual, is_na);
	RowMatrixXd p=cov.llt().matrixU();
	////cout << p << endl;

	vector<RowMatrixXd> MA_Phi_P;
	MA_Phi_P.reserve(ahead);
	for(int i0=0; i0<ahead; ++i0){
		RowMatrixXd temp_mat=p * ma_phi[i0];
		MA_Phi_P.push_back(temp_mat);

	}

	vector<RowMatrixXd> tbr;
	int num_dep=residual.cols();
	tbr.reserve(num_dep);
#pragma omp parallel for
	for(int i0=0; i0<num_dep;++i0){
		RowMatrixXd temp_mat(ahead, num_dep);
		for(int i=0;i<ahead;++i)
			temp_mat.row(i)=MA_Phi_P[i].row(i0);

		//cout <<temp_mat << endl;
		tbr.push_back(temp_mat);
	}

	return (tbr);
}


vector<RowMatrixXd>girf(RowMatrixXd &residual, VectorXi &is_na, RowMatrixXd &beta, int ahead, int lags){
	vector<RowMatrixXd> ma_phi=IRF_matrix(beta, ahead, lags);
	int num_dep=residual.cols();
	RowMatrixXd cov= vcov_residual(residual, is_na);


	vector<double> sigmas(num_dep);
	for(int i=0;i<num_dep;++i)
		sigmas[i]=1/sqrt(cov(i,i));

	vector<RowMatrixXd> MA_Phi_P(ahead);

	for(int i0=0; i0<ahead; ++i0){

		RowMatrixXd temp_mat=cov*ma_phi[i0];
		for (int r=0; r<temp_mat.rows(); ++r)
			for (int col=0; col<temp_mat.cols();++col)
				temp_mat(r,col)=temp_mat(r,col)*sigmas[r];

		MA_Phi_P[i0]=temp_mat;

	}

	vector<RowMatrixXd> tbr(num_dep);
//	tbr.reserve(num_dep);
	for(int i0=0; i0<num_dep;++i0){
		RowMatrixXd temp_mat(ahead, num_dep);
		for(int i=0;i<ahead;++i)
			temp_mat.row(i)=MA_Phi_P[i].row(i0);

		//cout <<temp_mat << endl;
		//cout <<"--"<<endl;
		tbr[i0]=temp_mat;
	}

	return (tbr);
}