#include "dual.hpp"
#include "quadrature_multipoly.hpp"
#include <stdexcept>
#include <vector>

constexpr int spatial_dim = 2;

template <typename T, int num_pts_1d> auto algoim_weight(std::vector<T> dv) {
  if (dv.size() != num_pts_1d * num_pts_1d) {
    throw std::runtime_error("incompatible design variable size");
  }
  int ndv = dv.size();

  std::vector<T> pts, wts;

  std::vector<T> data(num_pts_1d * num_pts_1d, 0.0);

  algoim::xarray<T, spatial_dim> phi(
      data.data(), algoim::uvector<int, spatial_dim>(num_pts_1d, num_pts_1d));

  algoim::bernstein::bernsteinInterpolate<spatial_dim>(
      [&](const algoim::uvector<T, spatial_dim> &xi) { // xi in [0, 1]
        return xi(0) + 0.9 * xi(1) - 1.5;
      },
      phi);

  for (int i = 0; i < ndv; i++) {
    data[i] += dv[i];
  }

  algoim::ImplicitPolyQuadrature<spatial_dim, T> ipquad(phi);
  ipquad.integrate(algoim::AutoMixed, num_pts_1d,
                   [&](const algoim::uvector<T, spatial_dim> &x, T w) {
                     if (algoim::bernstein::evalBernsteinPoly(phi, x) <= 0.0) {
                       for (int d = 0; d < spatial_dim; d++) {
                         pts.push_back(x(d));
                       }
                       wts.push_back(w);
                     }
                   });

  return std::make_tuple(pts, wts);
}

int main() {
  int constexpr num_pts_1d = 4;
  int j = 31; // quad index

  using T = double;
  using T2 = duals::dual<T>;

  double dh = 1e-6;

  int ndv = num_pts_1d * num_pts_1d;
  std::vector<T2> dv1(ndv, 0.0), dv2(ndv, 0.0);

  std::vector<T> w, grad_exact, grad_fd;

  for (int i = 0; i < ndv; i++) {
    dv1[i].dpart(1.0);
    dv2[i].rpart(dh);

    auto [pts1, w1] = algoim_weight<T2, num_pts_1d>(dv1);
    auto [pts2, w2] = algoim_weight<T2, num_pts_1d>(dv2);
    dv1[i].dpart(0.0);
    dv2[i].rpart(0.0);

    w.push_back(w1[j].rpart());
    grad_exact.push_back(w1[j].dpart());
    grad_fd.push_back((w2[j].rpart() - w1[j].rpart()) / dh);
  }

  for (int i = 0; i < ndv; i++) {
    T abserr = fabs(grad_exact[i] - grad_fd[i]);
    T relerr = abserr / fabs(grad_exact[i]);
    printf("w[%2d]:%15.5e, grad_exact[%2d]:%15.5e, grad_fd[%2d]:%15.5e, "
           "relerr:%15.5e\n",
           i, w[i], i, grad_exact[i], i, grad_fd[i], relerr);
  }
}
