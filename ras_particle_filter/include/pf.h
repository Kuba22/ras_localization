#include <cmath>
#include <armadillo>

using namespace arma;

class ParticleFilter {

	void init(vec bound, vec start_pose, mat& S, mat& R, mat& Q, double& Lambda_psi)
	{
		int M, part_bound;
		M = 1000;
		part_bound = 20;
		if (!start_pose.is_empty())
		{
			S = join_cols(repmat(start_pose, 1, M), (1.0 / M)*arma::ones<mat>(1, M));
		}
		else
		{
			S = arma::join_cols(arma::join_cols(arma::join_cols(arma::randu<mat>(1, M)*(bound(1) - bound(0) + 2 * part_bound) + bound(0) - part_bound, arma::randu<mat>(1, M)*(bound(3) - bound(2) + 2 * part_bound) + bound(2) - part_bound), arma::randu<mat>(1, M) * 2 * datum::pi - datum::pi), 1 * 1.0 / M*arma::ones<mat>(1, M));
		}
		vec R_diag = { 1e-2, 1e-2, 1e-2 };
		vec Q_diag = { 1e-1, 1e-1 };
		R = diagmat(R_diag);
		Q = diagmat(Q_diag);
		Lambda_psi = 0.0001;
	}

	mat observation_model(mat S, mat W, int j)
	{
		rowvec bearings = atan2(W(1, j - 1) - S.row(1), W(0, j - 1) - S.row(0)) - S.row(2) + datum::pi - floor((atan2(W(1, j - 1) - S.row(1), W(0, j - 1) - S.row(0)) - S.row(2) + datum::pi) / (2 * datum::pi)) * 2 * datum::pi - datum::pi;
		rowvec ranges = sqrt(arma::square((W(0, j - 1) - S.row(0))) + arma::square((W(1, j - 1) - S.row(1))));
		mat h = join_cols(ranges, bearings);
		return h;
	}

	mat predict(mat S, double v, double omega, mat R, double delta_t)
	{
		mat S_bar;
		rowvec _aux_rowvec_1;
		uword M;
		vec u;
		M = S.n_cols;
		u = arma::join_cols(arma::join_cols(arma::join_cols(v*arma::cos(S.row(2)), v*arma::sin(S.row(2))), omega*arma::ones<mat>(1, M)), arma::zeros<mat>(1, M))*delta_t;
		double __aux_rowvec_1[] = { 0, 0, 0 };
		_aux_rowvec_1 = rowvec(__aux_rowvec_1, 3, false);
		S_bar = S + u + {arma::trans(mvnrnd(_aux_rowvec_1, R, M)), arma::zeros<umat>(1, M)};
		return S_bar;
	}

	void associate(mat S_bar, mat z, mat W, double Lambda_psi, mat Q, rowvec& outlier, cube& Psi)
	{
		int M, N, c, i, k, m, n, psiSum;
		cube z_hat;
		cube z_m(2, M, N);
		double gaussian_multiplier, psi;
		mat D, Q_i;
		n = z.n_cols;
		N = W.n_cols;
		M = S_bar.n_cols;
		Psi = arma::zeros<cube>(1, n, M);
		outlier = arma::zeros<rowvec>(n);
		D = arma::zeros<mat>(M, N);
		cube v[2];
		v[0] = arma::zeros<cube>(n, M, N);
		v[1] = arma::zeros<cube>(n, M, N);
		z_hat = arma::zeros<cube>(2, M, N);
		gaussian_multiplier = (2 * datum::pi*sqrt(det(Q))) - 1;
		Q_i = Q - 1;
		for (k = 1; k <= N; k++)
		{
			//z_hat(span(0, z_hat.n_rows - 1), span(0, z_hat.n_cols - 1), k - 1) = observation_model(S_bar, W, k);
			z_hat.slice(k - 1) = observation_model(S_bar, W, k);
		}
		for (i = 1; i <= n; i++)
		{
			z_m.each_slice() = repmat(z.col(i - 1), 1, M);
			v(m2cpp::span<uvec>(0, v.n_rows - 1), i, m2cpp::span<uvec>(0, v.n_slices - 1), span::all) = double(arma::as_scalar(z_m(span(0, z_m.n_rows - 1), span(0, z_m.n_cols - 1), span(0, z_m.n_slices - 1)) - z_hat(span(0, z_hat.n_rows - 1), span(0, z_hat.n_cols - 1), span(0, z_hat.n_slices - 1))));
			v(2, i, m2cpp::span<uvec>(0, v.n_slices - 1), span::all) = v(2, i, m2cpp::span<uvec>(0, v.n_slices - 1), span::all) + datum::pi % 2 * datum::pi - datum::pi;
			for (k = 1; k <= N; k++)
			{
				D.col(k - 1) = arma::strans(arma::sum(arma:vectorize(reshape(v(m2cpp::span<uvec>(0, v.n_rows - 1), i, m2cpp::span<uvec>(0, v.n_slices - 1), k), { 2, M }) % (Q_i*reshape(v(m2cpp::span<uvec>(0, v.n_rows - 1), i, m2cpp::span<uvec>(0, v.n_slices - 1), k), { 2, M }))), 0));
			}
			psiSum = 0;
			D.min(c);
			for (m = 1; m <= M; m++)
			{
				psi = gaussian_multiplier*std::exp(-0.5*D(m - 1, c(m) - 1));
				Psi(0, i - 1, m - 1) = psi;
				psiSum = psiSum + psi;
			}
			outlier(i - 1) = psiSum*1.0 / M <= Lambda_psi;
		}
	}

	mat weight(mat S_bar, cube Psi, rowvec outlier)
	{
		mat p;
		p = prod(Psi(m2cpp::span<uvec>(0, Psi.n_rows - 1) - 1, outlier == 0 - 1, m2cpp::span<uvec>(0, Psi.n_slices - 1) - 1), 2);
		S_bar.row(3) = p / double(arma::as_scalar(arma::sum(arma:vectorize(p))));
		return S_bar;
	}

	mat systematic_resample(mat S_bar)
	{
		double r_0;
		int M, m;
		mat S, cdf;
		urowvec _aux_urowvec_1;
		uvec i;
		cdf = cumsum(S_bar.row(3));
		M = S_bar.row(0).n_cols;
		S.rows(arma::span(0, 2)).fill(arma::zeros<mat>(S_bar.rows(arma::span(0, 2)).n_rows, S_bar.rows(arma::span(0, 2)).n_cols));
		r_0 = double(arma::as_scalar(arma::randu<vec>(1) / M));
		for (m = 1; m <= M; m++)
		{
			i = find(cdf >= r_0, 1, "first") + 1;
			S(m2cpp::span<uvec>(0, 2), m2cpp::span<uvec>(m - 1, m - 1)) = S_bar(m2cpp::span<uvec>(0, 2), i);
			r_0 = r_0 + 1 * 1.0 / M;
		}
		_aux_urowvec_1 = { S_bar.row(3).n_rows, S_bar.row(3).n_cols };
		S.row(3).fill(1 * 1.0 / M*arma::ones<double>(_aux_urowvec_1));
		return S;
	}

	void mcl(mat S, mat R, mat Q, mat z, rowvec known_associations, double v, double omega, mat W, double Lambda_psi, rowvec Map_IDS, double delta_t, int t, mat& S, double& outliers)
	{
		cube Psi;
		int RESAMPLE_MODE, USE_KNOWN_ASSOCIATIONS, i;
		mat S_bar;
		rowvec outlier;
		uvec map_ids;
		S_bar = predict(S, v, omega, R, delta_t);
			associate(S_bar, z, W, Lambda_psi, Q, outlier, Psi);
		outliers = double(arma::as_scalar(arma::sum(outlier)));
		if (outliers)
		{
		}
		S_bar = weight(S_bar, Psi, outlier);
		RESAMPLE_MODE = 2;
		if (0 == RESAMPLE_MODE)
		{
			S = S_bar;
		}
		else if (2 == RESAMPLE_MODE)
		{
			S = systematic_resample(S_bar);
		}
	}
};