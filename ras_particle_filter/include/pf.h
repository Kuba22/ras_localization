#include <cmath>
#include <armadillo>
#include "PfHelper.h"

using namespace arma;

class ParticleFilter {
	default_random_engine gen;
	normal_distribution<double> normal_kv;
	normal_distribution<double> normal_kw;
	normal_distribution<double> normal_kd;
public:
	void init(vec bound, double part_bound, vec start_pose, mat& S, mat R, mat Q, int M, double particle_spread)
	{
		gen = default_random_engine(random_device()());
		normal_kd = normal_distribution<double>(0, sqrt(R(0, 0)) + 1e-9);
		normal_kv = normal_distribution<double>(0, sqrt(R(1, 1)) + 1e-9);
		normal_kw = normal_distribution<double>(0, sqrt(R(2, 2)) + 1e-9);
		if (!start_pose.is_empty() && particle_spread>0)
		{
			normal_distribution<double> normal(0, particle_spread + 1e-9);
			std::uniform_real_distribution<double> unif(0.0, 2*datum::pi);
			S = mat(4, M);
			double iM = 1.0 / M;
			for (int m = 0; m<M; m++)
			{
				S(0, m) = start_pose(0) + normal(gen);
				S(1, m) = start_pose(1) + normal(gen);
				S(2, m) = unif(gen);
				S(3, m) = iM;
			}
		}
		else if (!start_pose.is_empty())
		{
			S = join_cols(repmat(start_pose, 1, M), (1.0 / M)*arma::ones<mat>(1, M));
		}
		else
		{
			S = arma::join_cols(arma::join_cols(arma::join_cols(arma::randu<mat>(1, M)*(bound(1) - bound(0) + 2 * part_bound) + bound(0) - part_bound, arma::randu<mat>(1, M)*(bound(3) - bound(2) + 2 * part_bound) + bound(2) - part_bound), arma::randu<mat>(1, M) * 2 * datum::pi - datum::pi), 1.0 / M*arma::ones<mat>(1, M));
		}
	}

	mat observation_model(mat S, mat W, double phi)
	{
		int M = S.n_cols;
		mat h(2, M);
		for (int m = 0; m < M; m++) {
			h.col(m) = vec({ getRange(W, S.col(m), phi), phi });
		}
		return h;
	}

	mat predict(mat S, double v, double omega, double delta_t)
	{
		int M = S.n_cols;

		mat dS(4, M);
		for (int m = 0; m < M; m++) {
			dS.col(m) = vec({
				cos(S(2, m))*(v*delta_t)*(v*delta_t)*normal_kd(gen),
				sin(S(2, m))*(v*delta_t)*(v*delta_t)*normal_kd(gen),
				(v*delta_t)*(v*delta_t)*normal_kv(gen) + (omega*delta_t)*(omega*delta_t)*normal_kw(gen), 0 });
		}
		mat u = join_cols(join_cols(join_cols(v*arma::cos(S.row(2)), v*arma::sin(S.row(2))), omega*ones<mat>(1, M)), zeros<mat>(1, M))*delta_t;
		mat S_bar = S + u + dS;
		return S_bar;
	}

	void associate(mat S_bar, mat z, mat W, double Lambda_psi, mat Q, rowvec& outlier, cube& Psi)
	{
		int n = z.n_cols;
		int N = z.n_cols;
		int M = S_bar.n_cols;
		Psi = arma::zeros<cube>(1, n, M);
		outlier = arma::zeros<rowvec>(n);
		cube z_hat = arma::zeros<cube>(2, M, N);
		double gaussian_multiplier = 1 / (2 * datum::pi*sqrt(det(Q)));
		mat Q_i = Q.i();
		for (int k = 0; k < N; k++)
		{
			z_hat.slice(k) = observation_model(S_bar, W, z(1, k));
		}
		for (int i = 0; i < n; i++)
		{
			double psiSum = 0;
			for (int m = 0; m < M; m++)
			{
				vec v = z.col(i) - z_hat.slice(i).col(m);
				v(1) = mod(v(1) + datum::pi, 2 * datum::pi) - datum::pi;
				double D = as_scalar(v.t()*Q_i*v);
				double psi = gaussian_multiplier*std::exp(-0.5*D);
				Psi(0, i, m) = psi;
				psiSum = psiSum + psi;
			}
			outlier(i) = psiSum*1.0 / M <= Lambda_psi;
		}
	}

	mat weight(mat S_bar, cube Psi, rowvec outlier)
	{
		int n = outlier.size();
		int M = S_bar.n_cols;

		mat psi = reshape(Psi, n, M, 1);
		rowvec p = prod(psi.rows(find(outlier == false)), 0);
		S_bar.row(3) = p / sum(p);
		return S_bar;
	}

	mat resample(mat S_bar)
	{
		rowvec cdf = cumsum(S_bar.row(3));
		int M = S_bar.n_cols;
		mat S(4, M);
		S.rows(arma::span(0, 2)).fill(0);
		double r_0 = double(arma::as_scalar(arma::randu<vec>(1) / M));
		for (int m = 0; m < M; m++)
		{
			int i = as_scalar(find(cdf >= r_0, 1));
			S(0, m) = S_bar(0, i);
			S(1, m) = S_bar(1, i);
			S(2, m) = S_bar(2, i);
			r_0 = r_0 + 1.0 / M;
		}
		S.row(3) = 1.0 / M * ones(1, M);
		return S;
	}
};