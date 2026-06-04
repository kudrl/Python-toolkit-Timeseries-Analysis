"""Скользящее окно и 3D-сканирование (window × lag × position) для анализа во времени."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from .statistics import lag_quality, select_best_median_worst


def analyze_sliding_windows(
    data: pd.DataFrame,
    variant: str,
    window_size: int,
    stride: int,
    *,
    compute_variant_func,
    is_pvalue: bool,
    lag: int = 1,
    pairs: list[tuple[int, int]] | None = None,
    start_min: int | None = None,
    start_max: int | None = None,
    max_windows: int = 400,
    return_matrices: bool = False,
    n_jobs: int | None = None,
    parallel_backend: str | None = None,
) -> dict:
    """Скользящее окно заданного размера.

    Возвращает структуру для HTML-отчёта:
      {
        "best_window": {"start": int, "end": int, "metric": float, "matrix": ndarray},
        "curve": {"x": [start_idx...], "y": [metric...]},
        "ticks": [{"start":.., "end":.., "metric":.., "matrix":..}, ...],
        "extremes": {"best": idx, "median": idx, "worst": idx}
      }
    """
    if data is None or data.empty:
        return {}

    n = len(data)
    w = int(max(2, min(window_size, n)))
    s = int(max(1, stride))

    best = {"start": 0, "end": w, "metric": float("-inf"), "matrix": None}

    st0 = int(max(0, start_min)) if start_min is not None else 0
    st1 = int(min(n - w, start_max)) if start_max is not None else (n - w)
    st1 = int(max(st0, st1))

    max_windows = int(max(1, max_windows))
    starts = list(range(st0, st1 + 1, s))
    if len(starts) > max_windows:
        idx = np.linspace(0, len(starts) - 1, max_windows).round().astype(int)
        starts = [starts[i] for i in idx]

    try:
        nj = int(n_jobs) if n_jobs is not None else 1
    except (TypeError, ValueError, OverflowError):
        nj = 1
    nj = int(max(1, nj))

    def _compute_one(start: int) -> tuple[int, int, float, np.ndarray | None]:
        end = start + w
        if end > n:
            return int(start), int(end), float("nan"), None
        chunk = data.iloc[start:end]
        try:
            mat = compute_variant_func(chunk, variant, lag=int(max(1, lag)), pairs=pairs)
            score = lag_quality(variant, mat, is_pvalue)
            return (
                int(start),
                int(end),
                float(score) if np.isfinite(score) else float("nan"),
                (mat if return_matrices else None),
            )
        except Exception as ex:
            logging.error("[SlidingWindow] %s win=%d start=%d: %s", variant, w, start, ex)
            return int(start), int(end), float("nan"), None

    if nj == 1 or len(starts) <= 1:
        computed = [_compute_one(st) for st in starts]
    else:
        try:
            from joblib import Parallel, delayed

            backend = str(parallel_backend or "loky")
            computed = Parallel(n_jobs=nj, backend=backend)(
                delayed(_compute_one)(st) for st in starts
            )
        except ImportError:
            computed = [_compute_one(st) for st in starts]

    xs: list[int] = []
    ys: list[float] = []
    ticks: list[dict] = []
    for start, end, score_f, mat in computed:
        xs.append(int(start))
        ys.append(float(score_f) if np.isfinite(score_f) else float("nan"))
        ticks.append(
            {
                "start": int(start),
                "end": int(end),
                "metric": float(score_f) if np.isfinite(score_f) else float("nan"),
                "matrix": mat if return_matrices else None,
            }
        )
        if np.isfinite(score_f) and float(score_f) > float(best["metric"]):
            best = {"start": int(start), "end": int(end), "metric": float(score_f), "matrix": mat}

    return {
        "best_window": best,
        "curve": {"x": xs, "y": ys},
        "ticks": ticks,
        "extremes": select_best_median_worst(ticks, key="metric"),
    }


def analyze_window_lag_cube(
    data: pd.DataFrame,
    variant: str,
    *,
    window_sizes: list[int],
    lag_grid: list[int],
    stride: int,
    compute_variant_func,
    is_pvalue: bool,
    eval_limit: int = 120,
    matrix_limit: int = 60,
    max_windows_per_size: int = 80,
) -> dict:
    """3D-скан window_size × lag × start_pos для отчета.

    Возвращает компактный payload, совместимый с HTML-рендером:
    ``points`` для scatter3d, ``gallery`` и ``matrices`` для просмотра теплокарт.
    """
    if data is None or data.empty:
        return {}

    n = int(len(data))
    sizes = sorted({int(max(2, min(w, n))) for w in (window_sizes or []) if int(w) > 1})
    lags = sorted({int(max(1, lag)) for lag in (lag_grid or [])})
    if not sizes or not lags:
        return {}

    stride = int(max(1, stride))
    eval_limit = int(max(1, eval_limit))
    matrix_limit = int(max(0, matrix_limit))
    max_windows_per_size = int(max(1, max_windows_per_size))

    jobs: list[tuple[int, int, int, int]] = []
    for w in sizes:
        starts = list(range(0, max(0, n - w) + 1, stride))
        if len(starts) > max_windows_per_size:
            idx = np.linspace(0, len(starts) - 1, max_windows_per_size).round().astype(int)
            starts = [starts[i] for i in idx]
        for lag in lags:
            for start in starts:
                jobs.append((int(w), int(lag), int(start), int(start + w)))

    if len(jobs) > eval_limit:
        idx = np.linspace(0, len(jobs) - 1, eval_limit).round().astype(int)
        jobs = [jobs[i] for i in idx]

    points: list[dict] = []
    matrices_added = 0
    for idx, (w, lag, start, end) in enumerate(jobs):
        pid = f"w{w}_l{lag}_s{start}"
        try:
            chunk = data.iloc[start:end]
            mat = compute_variant_func(chunk, variant, lag=int(lag))
            score = lag_quality(variant, mat, is_pvalue)
            score_f = float(score) if np.isfinite(score) else float("nan")
            matrix = mat if matrices_added < matrix_limit else None
            if matrix is not None:
                matrices_added += 1
        except Exception as ex:
            logging.error("[WindowLagCube] %s w=%d lag=%d start=%d: %s", variant, w, lag, start, ex)
            score_f = float("nan")
            matrix = None
        points.append(
            {
                "id": pid,
                "window_size": int(w),
                "lag": int(lag),
                "start": int(start),
                "end": int(end),
                "metric": score_f,
                "matrix": matrix,
            }
        )

    idxs = select_best_median_worst(points, key="metric")
    id_extremes = {
        key: (points[val]["id"] if isinstance(val, int) and 0 <= val < len(points) else None)
        for key, val in idxs.items()
    }
    for tag, pid in id_extremes.items():
        if pid is None:
            continue
        for point in points:
            if point.get("id") == pid:
                point["tag"] = tag
                break

    gallery_ids = [pid for pid in [id_extremes.get("best"), id_extremes.get("median"), id_extremes.get("worst")] if pid]
    gallery = [dict(p) for p in points if p.get("id") in set(gallery_ids)]
    selectable_ids = [str(p["id"]) for p in points if p.get("matrix") is not None]

    return {
        "matrix_mode": "limited",
        "matrix_limit": int(matrix_limit),
        "eval_limit": int(eval_limit),
        "combos": [{"window_size": w, "lag": lag} for w in sizes for lag in lags],
        "window_sizes": sizes,
        "lag_grid": lags,
        "points": points,
        "gallery": gallery,
        "selectable_ids": selectable_ids,
        "extremes": id_extremes,
    }
