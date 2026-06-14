use extendr_api::prelude::*;

use runorm::{
    estimate_overdispersion, normalize_csr, CsrCounts, NormParams, PfTarget, ShiftedClrMatrix,
};
use rupca::{
    pca_scanpy_sparse_csr, pca_shifted_clr_sparse_csr, CsrMatrix, ScanpyPcaParams,
    ScanpyPcaResult, ShiftedClrCsrMatrix,
};

fn norm_err(e: runorm::NormError) -> Error {
    Error::Other(e.to_string())
}

fn pca_err(e: rupca::RuPcaError) -> Error {
    Error::Other(e.to_string())
}

fn build_counts(
    data: Doubles,
    indices: Integers,
    indptr: Integers,
    n_cells: i32,
    n_features: i32,
) -> Result<CsrCounts> {
    let data = data.iter().map(|x| x.inner()).collect::<Vec<f64>>();
    let indices = indices.iter().map(|x| x.inner() as usize).collect::<Vec<usize>>();
    let indptr = indptr.iter().map(|x| x.inner() as usize).collect::<Vec<usize>>();
    CsrCounts::new(n_cells as usize, n_features as usize, data, indices, indptr).map_err(norm_err)
}

fn build_sparse(
    data: Doubles,
    indices: Integers,
    indptr: Integers,
    n_cells: i32,
    n_features: i32,
) -> CsrMatrix {
    CsrMatrix {
        n_rows: n_cells as usize,
        n_cols: n_features as usize,
        data: data.iter().map(|x| x.inner()).collect::<Vec<f64>>(),
        indices: indices.iter().map(|x| x.inner() as usize).collect::<Vec<usize>>(),
        indptr: indptr.iter().map(|x| x.inner() as usize).collect::<Vec<usize>>(),
    }
}

fn parse_target(target: &str, fixed: Robj, alpha: Robj) -> Result<PfTarget> {
    if !alpha.is_null() {
        let a = alpha.as_real().ok_or_else(|| Error::Other("alpha must be numeric".into()))?;
        return Ok(PfTarget::Alpha(a));
    }
    match target {
        "mean" => Ok(PfTarget::MeanDepth),
        "median" => Ok(PfTarget::MedianDepth),
        "auto" => Ok(PfTarget::EstimateAlpha),
        "fixed" => {
            let k = fixed
                .as_real()
                .ok_or_else(|| Error::Other("target='fixed' requires a numeric K".into()))?;
            Ok(PfTarget::Fixed(k))
        }
        other => Err(Error::Other(format!(
            "unknown target '{other}' (expected mean, median, auto, or fixed)"
        ))),
    }
}

fn int_vec(v: &[usize]) -> Vec<i32> {
    v.iter().map(|&x| x as i32).collect()
}

fn pca_params(n_components: i32, ncv: Robj, maxiter: Robj, seed: f64, tol: f64) -> ScanpyPcaParams {
    ScanpyPcaParams {
        n_components: n_components as usize,
        tol,
        ncv: if ncv.is_null() { None } else { ncv.as_integer().map(|x| x as usize) },
        maxiter: if maxiter.is_null() { None } else { maxiter.as_integer().map(|x| x as usize) },
        seed: seed as u64,
    }
}

fn pca_result_list(r: ScanpyPcaResult) -> List {
    list!(
        scores = r.scores,
        components = r.components,
        mean = r.mean,
        explained_variance = r.explained_variance,
        explained_variance_ratio = r.explained_variance_ratio,
        singular_values = r.singular_values,
        noise_variance = r.noise_variance,
        n_samples = r.n_samples as i32,
        n_features = r.n_features as i32,
        n_components = r.n_components as i32
    )
}

fn pca_result_list_with_norm(r: ScanpyPcaResult, k_value: f64, alpha: Option<f64>) -> List {
    list!(
        scores = r.scores,
        components = r.components,
        mean = r.mean,
        explained_variance = r.explained_variance,
        explained_variance_ratio = r.explained_variance_ratio,
        singular_values = r.singular_values,
        noise_variance = r.noise_variance,
        n_samples = r.n_samples as i32,
        n_features = r.n_features as i32,
        n_components = r.n_components as i32,
        k = k_value,
        alpha = alpha
    )
}

/// Estimate overdispersion alpha and K from a features-by-cells dgCMatrix slot representation.
/// The R layer passes dgC slots directly; they are interpreted as CSR over cells.
#[extendr]
fn scclr_overdispersion(
    data: Doubles,
    indices: Integers,
    indptr: Integers,
    n_cells: i32,
    n_features: i32,
) -> Result<List> {
    let counts = build_counts(data, indices, indptr, n_cells, n_features)?;
    let od = estimate_overdispersion(&counts).map_err(norm_err)?;
    Ok(list!(alpha = od.alpha, mean_depth = od.mean_depth, k = od.k))
}

/// PFlogPF normalization. Returns sparse shifted-log CSR slots plus the per-cell center vector.
#[extendr]
fn scclr_normalize(
    data: Doubles,
    indices: Integers,
    indptr: Integers,
    n_cells: i32,
    n_features: i32,
    target: &str,
    fixed: Robj,
    alpha: Robj,
    log1p: bool,
    center: bool,
) -> Result<List> {
    let counts = build_counts(data, indices, indptr, n_cells, n_features)?;
    let params = NormParams {
        target: parse_target(target, fixed, alpha)?,
        log1p,
        center,
    };
    let (m, report) = normalize_csr(&counts, &params).map_err(norm_err)?;
    Ok(list!(
        data = m.data,
        indices = int_vec(&m.indices),
        indptr = int_vec(&m.indptr),
        row_center = m.row_center,
        k = report.k,
        alpha = report.alpha
    ))
}

/// Sparse PCA. If row_center is non-NULL, PCA is run on sparse - row_center without densifying.
#[extendr]
fn scclr_pca(
    data: Doubles,
    indices: Integers,
    indptr: Integers,
    n_cells: i32,
    n_features: i32,
    row_center: Robj,
    n_components: i32,
    ncv: Robj,
    maxiter: Robj,
    seed: f64,
    tol: f64,
) -> Result<List> {
    let sparse = build_sparse(data, indices, indptr, n_cells, n_features);
    let params = pca_params(n_components, ncv, maxiter, seed, tol);
    let result = if row_center.is_null() {
        pca_scanpy_sparse_csr(&sparse, params).map_err(pca_err)?
    } else {
        let center = row_center
            .as_real_vector()
            .ok_or_else(|| Error::Other("row_center must be numeric or NULL".into()))?;
        let shifted = ShiftedClrCsrMatrix { sparse, row_center: center };
        pca_shifted_clr_sparse_csr(&shifted, params).map_err(pca_err)?
    };
    Ok(pca_result_list(result))
}

/// One-shot raw counts -> PFlogPF -> sparse PCA.
#[extendr]
fn scclr_normalize_pca(
    data: Doubles,
    indices: Integers,
    indptr: Integers,
    n_cells: i32,
    n_features: i32,
    n_components: i32,
    target: &str,
    fixed: Robj,
    alpha: Robj,
    ncv: Robj,
    maxiter: Robj,
    seed: f64,
    tol: f64,
) -> Result<List> {
    let counts = build_counts(data, indices, indptr, n_cells, n_features)?;
    let params = NormParams {
        target: parse_target(target, fixed, alpha)?,
        log1p: true,
        center: true,
    };
    let (m, report) = normalize_csr(&counts, &params).map_err(norm_err)?;
    let ShiftedClrMatrix { n_rows, n_cols, data, indices, indptr, row_center } = m;
    let shifted = ShiftedClrCsrMatrix {
        sparse: CsrMatrix { n_rows, n_cols, data, indices, indptr },
        row_center,
    };
    let result = pca_shifted_clr_sparse_csr(
        &shifted,
        pca_params(n_components, ncv, maxiter, seed, tol),
    )
    .map_err(pca_err)?;
    Ok(pca_result_list_with_norm(result, report.k, report.alpha))
}

extendr_module! {
    mod scclrR;
    fn scclr_overdispersion;
    fn scclr_normalize;
    fn scclr_pca;
    fn scclr_normalize_pca;
}
