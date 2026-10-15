"""Small OLS helper replacing MATLAB fitlm (coefficients, p-values, R^2, prediction bands)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass
class OLSResult:
    names: list
    beta: np.ndarray
    se: np.ndarray
    pvalues: np.ndarray
    r2: float
    sigma2: float
    XtX_inv: np.ndarray
    dof: int

    def predict(self, X, interval: str | None = "observation", alpha: float = 0.05):
        """Predicted mean and (optional) prediction interval, like fitlm predict(...,'Prediction','observation')."""
        X = np.atleast_2d(np.asarray(X, float))
        if X.shape[1] == len(self.beta) - 1:
            X = np.column_stack([np.ones(len(X)), X])
        y = X @ self.beta
        if interval is None:
            return y
        v = np.einsum("ij,jk,ik->i", X, self.XtX_inv, X) * self.sigma2
        if interval == "observation":
            v = v + self.sigma2
        h = stats.t.ppf(1 - alpha / 2, self.dof) * np.sqrt(v)
        return y, np.column_stack([y - h, y + h])

    def summary(self) -> str:
        rows = [f"{'term':>12} {'estimate':>12} {'SE':>12} {'p':>10}"]
        for n, b, s, p in zip(self.names, self.beta, self.se, self.pvalues):
            rows.append(f"{n:>12} {b:12.4g} {s:12.4g} {p:10.3g}")
        rows.append(f"R^2 = {self.r2:.4f}, dof = {self.dof}")
        return "\n".join(rows)


def ols(y, X, names=None, intercept=True) -> OLSResult:
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    names = list(names or [f"x{i + 1}" for i in range(X.shape[1])])
    if intercept:
        X = np.column_stack([np.ones(len(y)), X])
        names = ["(Intercept)"] + names
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    sigma2 = resid @ resid / dof if dof > 0 else np.nan
    XtX_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(XtX_inv) * sigma2)
    with np.errstate(divide="ignore", invalid="ignore"):
        tvals = beta / se
    p = 2 * stats.t.sf(np.abs(tvals), dof) if dof > 0 else np.full_like(beta, np.nan)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - resid @ resid / ss_tot if ss_tot > 0 else np.nan
    return OLSResult(names, beta, se, p, r2, sigma2, XtX_inv, dof)
