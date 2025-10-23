"""
Pedestrian Risk Plotter
-----------------------
Plots pedestrian risk vs. impact speed for:
  - Tefft (2013) death risk (U.S. crash-reconstruction logistic model)
  - Monfort & Mueller (2025) injury risk (fatal / MAIS 3+F / MAIS 2+F),
    with hood leading-edge (HLE) height, age, sex, and speed (km/h).

Notes
- Tefft model applies to pedestrians >= 15 years; this script guards and can fallback to an adult-average speed-only curve.
- Mueller model uses speeds internally in km/h; x-axis can be shown in mph or km/h.
"""

import argparse
import json
from typing import Union, Sequence, Optional, Dict, Any, Tuple, List

import numpy as np
import matplotlib.pyplot as plt
from math import atan2, degrees
import matplotlib.patheffects as pe

ArrayLike = Union[float, int, np.ndarray]

# -------------------- 
# Dark style (similar to plotly_dark) 
# --------------------
PLOTLY_DARK_RC = {
    "figure.facecolor": "#111111",
    "axes.facecolor":   "#111111",
    "axes.edgecolor":   "#444444",
    "axes.labelcolor":  "#E5E5E5",
    "xtick.color":      "#CCCCCC",
    "ytick.color":      "#CCCCCC",
    "text.color":       "#E5E5E5",
    "grid.color":       "#333333",
    "grid.alpha":       0.2,
    "axes.grid":        True,
    "legend.frameon":   False,
}

# -------------------- 
# Helpers 
# --------------------
def _x_at_y(xs: np.ndarray, ys: np.ndarray, yq: float):
    """Linear interpolation: return x where y crosses yq (if within range)."""
    if not (ys.min() <= yq <= ys.max()):
        return None
    idx = np.searchsorted(ys, yq)
    if idx == 0 or idx >= len(xs):
        return None
    x0, x1 = xs[idx-1], xs[idx]
    y0, y1 = ys[idx-1], ys[idx]
    if y1 == y0:
        return float(x0)
    return float(x0 + (yq - y0) * (x1 - x0) / (y1 - y0))

def _angle_at_x(xs: np.ndarray, ys: np.ndarray, xq: float):
    """Angle (deg) of curve near x==xq using local secant."""
    i = np.searchsorted(xs, xq)
    i0 = max(1, min(i, len(xs) - 1))  # clamp
    x0, x1 = xs[i0 - 1], xs[i0]
    y0, y1 = ys[i0 - 1], ys[i0]
    return degrees(atan2(y1 - y0, x1 - x0))

def _resolve_rotation(rotation_opt: Union[str, float, None], xs, ys, xq):
    """'auto' -> curve angle; float/int -> fixed degrees; None/other -> 0."""
    if rotation_opt is None:
        return 0.0
    if isinstance(rotation_opt, (int, float)):
        return float(rotation_opt)
    if isinstance(rotation_opt, str) and rotation_opt.lower() == "auto":
        return _angle_at_x(xs, ys, xq)
    return 0.0

# -------------------- 
# Tefft (2013) death-risk model 
# --------------------
def _bmi_from_imp(h_in: ArrayLike, w_lb: ArrayLike):
    h_m = np.asarray(h_in, dtype=float) * 0.0254
    w_kg = np.asarray(w_lb, dtype=float) * 0.45359237
    return w_kg / (h_m ** 2)

def ped_death_probability_tefft(
                speed_mph: ArrayLike,
                age_years: ArrayLike,
                height_in: ArrayLike,
                weight_lb: ArrayLike,
                bmi: ArrayLike = None,
                vehicle_is_truck: ArrayLike = 0):
    v = np.asarray(speed_mph, dtype=float)
    a = np.asarray(age_years, dtype=float)
    h = np.asarray(height_in, dtype=float)
    w = np.asarray(weight_lb, dtype=float)
    t = np.asarray(vehicle_is_truck, dtype=float)
    BMI = _bmi_from_imp(h, w) if bmi is None else np.asarray(bmi, dtype=float)
    BMI_above25 = np.maximum(BMI - 25.0, 0.0)
    w_term = (w / 100.0) ** (-0.5)
    z = (
        196.8
        + 0.202 * v
        + 0.00712 * (a / 10.0) ** 3
        - 1.174 * h
        - 90.5 * w_term
        - 2.247 * BMI
        + 1.557 * BMI_above25
        + 0.543 * t
    )
    return 1.0 / (1.0 + np.exp(-z))

def ped_death_probability_tefft_guarded(
                speed_mph: ArrayLike,
                age_years: ArrayLike,
                height_in: ArrayLike,
                weight_lb: ArrayLike,
                bmi: ArrayLike = None,
                vehicle_is_truck: ArrayLike = 0,
                fallback: str = "speed_only"):
    """
    Guard for <15 y/o (Tefft not valid). fallback='speed_only' uses adult-average
    speed-only curve
    """
    age_arr = np.asarray(age_years, dtype=float)
    if np.isscalar(age_arr):
        under = (age_arr < 15)
    else:
        under = np.any(age_arr < 15)

    if under:
        if fallback == "speed_only":
            alpha = -5.406504111506091
            beta  =  0.1331651258991648
            v = np.asarray(speed_mph, dtype=float)
            z = alpha + beta * v
            return 1.0 / (1.0 + np.exp(-z))
        raise ValueError("Tefft (2013) is not valid for age < 15; provide age >= 15.")
    return ped_death_probability_tefft(speed_mph, age_years, height_in, weight_lb, bmi, vehicle_is_truck)

# -------------------- 
# Monfort & Mueller (2025) injury-risk model 
# --------------------
# Mean-centering constants and coefficients (Table 4)
_MUELLER_MEANS = dict(speed_kmh=42.79, hle_cm=85.25, sex_male=0.59, age_years=42.42)
_MUELLER_COEFS = {
    "mais2p":  dict(b0= 0.686, bv=0.083, bh=0.034, bs= 0.047, ba=0.023, bvh=0.001),
    "mais3p":  dict(b0=-0.460, bv=0.091, bh=0.039, bs=-0.119, ba=0.029, bvh=0.002),
    "fatal":   dict(b0=-3.280, bv=0.123, bh=0.047, bs= 0.319, ba=0.030, bvh=0.000),
}

def mueller_injury_probability(
                    speed_kmh: ArrayLike,
                    hle_cm: ArrayLike,
                    sex_male: ArrayLike,
                    age_years: ArrayLike,
                    severity: str = "fatal"):
    if severity not in _MUELLER_COEFS:
        raise ValueError("severity must be one of: 'fatal', 'mais3p', 'mais2p'")
    c = _MUELLER_COEFS[severity]
    m = _MUELLER_MEANS
    v = np.asarray(speed_kmh, dtype=float)
    X1 = v - m["speed_kmh"]
    X2 = float(hle_cm) - m["hle_cm"]
    X3 = float(sex_male) - m["sex_male"]
    X4 = float(age_years) - m["age_years"]
    eta = c["b0"] + c["bv"]*X1 + c["bh"]*X2 + c["bs"]*X3 + c["ba"]*X4 + c["bvh"]*(X1*X2)
    return 1.0 / (1.0 + np.exp(-eta))

# ----------------
# vehicle HLE assumptions 
# ----------------
def vehicle_hle_guess(
                vehicle_type: str,
                *, 
                suv_as: str = "pickup", 
                van_as: str = "pickup"):
    """
    Guess hood leading-edge (HLE) height in cm for a vehicle type.

    Defaults (per Mueller/Monfort medians):
      - Passenger car ≈ 75 cm
      - Pickup truck  ≈ 109 cm
    SUVs and vans map to 'pickup' or 'car' per the suv_as/van_as options.
    """
    s = str(vehicle_type).strip().lower()

    car_aliases = {"car", "sedan", "hatchback", "coupe", "wagon", "passenger car", "saloon"}
    pickup_aliases = {"pickup", "pickup truck", "truck", "pick-up", "ute"}
    suv_aliases = {"suv", "sport utility vehicle", "crossover", "cuv"}
    van_aliases = {"van", "minivan", "mpv", "people carrier"}

    if s in car_aliases:
        return 75.0
    if s in pickup_aliases:
        return 109.0
    if s in suv_aliases:
        return 109.0 if suv_as == "pickup" else 75.0
    if s in van_aliases:
        return 109.0 if van_as == "pickup" else 75.0

    raise ValueError(f"Unknown vehicle_type '{vehicle_type}'. "
                     "Use car/sedan/etc., pickup/truck, suv/crossover, or van/minivan.")

# -------------------- 
# Labels 
# --------------------
def _label_tefft(age, h_in, w_lb, truck_flag, pid=None, include_id=False):
    veh = "Light truck/SUV" if int(truck_flag) == 1 else "Car"
    base = f"{int(age)}y, {int(h_in)}\" {int(w_lb)} lb, {veh}"
    return (f"[{pid}] Tefft — " + base) if (include_id and pid is not None) else ("Tefft — " + base)

def _label_mueller(age, hle_cm, sex_male, severity, pid=None, include_id=False):
    sex = "Male" if int(sex_male) == 1 else "Female"
    sev = {"fatal":"Fatal", "mais3p":"MAIS 3+F", "mais2p":"MAIS 2+F"}[severity]
    base = f"{int(age)}y, HLE {hle_cm:.0f} cm, {sex}, {sev}"
    return (f"[{pid}] Mueller — " + base) if (include_id and pid is not None) else ("Mueller — " + base)

# ----------------------------------------
# plotting function 
# ----------------------------------------
def plot_ped_risk(
    *,
    model: str = "tefft", # 'tefft' | 'mueller' | 'both'
    units: str = "mph", # 'mph' or 'kmh' (display units)
    speed_min: float = 5,
    speed_max: float = 80,
    speed_step: float = 5,
    # Tefft inputs (single or multi)
    tefft_profile: Optional[Dict[str, Any]] = None,
    tefft_profiles: Optional[Sequence[Dict[str, Any]]] = None,
    # Mueller inputs (single or multi)
    mueller_profile: Optional[Dict[str, Any]] = None,
    mueller_profiles: Optional[Sequence[Dict[str, Any]]] = None,
    mueller_severity: str = "fatal",
    # markers
    pct_markers: Sequence[float] = (10, 25, 50, 75),
    show_marker_labels: bool = False,
    marker_size: float = 48,
    marker_kwargs: Optional[Dict[str, Any]] = None,
    marker_label_rotation: Union[str, float, None] = "auto",
    marker_label_offset: Tuple[int, int] = (6, 4),
    use_text_stroke: bool = True,
    # style/save
    title: Optional[str] = None,
    out: Optional[str] = "ped_risk_plot.png",
    show: bool = True,
    ax: Optional[plt.Axes] = None,
    # color pairing
    match_colors: bool = True,
    match_by: str = "id", # 'id' or 'position'
    id_key: str = "profile_id",
    include_id_in_label: bool = False
    ):
    """
    Plot pedestrian risk vs. speed for Tefft (2013), Mueller (2025), or both.

    TEFFT profile dict keys:
        age_years, height_in, weight_lb, (bmi), (vehicle_is_truck), (profile_id)

    MUELLER profile dict keys:
        age_years, hle_cm OR vehicle_type, sex_male, (severity -> 'fatal'|'mais3p'|'mais2p'), (profile_id)
    """
    # Speed array in display units
    x = np.arange(speed_min, speed_max + speed_step, speed_step)
    # Mueller needs km/h internally:
    x_kmh = x * 1.609344 if units.lower() == "mph" else x
    x_label = "Impact speed (mph)" if units.lower() == "mph" else "Impact speed (km/h)"

    # Normalize inputs into lists
    t_list = tefft_profiles if tefft_profiles is not None else ([tefft_profile] if tefft_profile else [])
    m_list = mueller_profiles if mueller_profiles is not None else ([mueller_profile] if mueller_profile else [])
    if model == "tefft":
        m_list = []
    elif model == "mueller":
        t_list = []
    elif model != "both":
        raise ValueError("model must be 'tefft', 'mueller', or 'both'")

    q_levels = np.array(pct_markers, dtype=float) / 100.0
    if marker_kwargs is None:
        marker_kwargs = dict(edgecolors="white", linewidths=0.8, zorder=5)

    # ---------- 
    # Color pairing logic 
    # ----------
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', None) or [
        "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
        "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"
    ]
    color_map: Dict[Any, str] = {}
    next_color_idx = 0
    def _color_for_id(pid):
        nonlocal next_color_idx
        if pid not in color_map:
            color_map[pid] = color_cycle[next_color_idx % len(color_cycle)]
            next_color_idx += 1
        return color_map[pid]

    with plt.rc_context(PLOTLY_DARK_RC):
        created_fig = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(9, 5.5))
            created_fig = True
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        # -------- 
        # Plot TEFFT (solid) 
        # --------
        tefft_lines: List[Tuple[plt.Line2D, Any]] = []
        for i, p in enumerate(t_list):
            if p is None:
                continue
            age = p["age_years"]; h = p["height_in"]; w = p["weight_lb"]
            b = p.get("bmi", None)
            tflag = p.get("vehicle_is_truck", 0)
            pid = p.get(id_key, None)

            y = ped_death_probability_tefft_guarded(x, age, h, w, b, tflag, fallback="speed_only")

            # Color selection
            if match_colors:
                if match_by == "id" and pid is not None:
                    color = _color_for_id(("pair", pid))
                elif match_by == "position":
                    color = color_cycle[i % len(color_cycle)]
                else:
                    color = None
            else:
                color = None

            (line,) = ax.plot(x, y, linewidth=2.6, linestyle='-',
                              label=_label_tefft(age, h, w, tflag, pid, include_id_in_label),
                              color=color)
            tefft_lines.append((line, pid))

            # Dots and labels
            for q in q_levels:
                xq = _x_at_y(x, y, q)
                if xq is None:
                    continue
                ax.scatter([xq], [q], s=marker_size, color=line.get_color(), **marker_kwargs)
                if show_marker_labels:
                    angle = _resolve_rotation(marker_label_rotation, x, y, xq)
                    ann = ax.annotate(
                        f"  {int(round(q*100))}% @ {xq:.1f} {units}",
                        xy=(xq, q), xytext=marker_label_offset, textcoords="offset points",
                        va="bottom", ha="left", fontsize=9,
                        rotation=angle, rotation_mode="anchor",
                    )
                    if use_text_stroke:
                        ann.set_path_effects([pe.withStroke(linewidth=3, foreground="black", alpha=0.6)])

        # Helper to find paired Tefft color
        def _paired_tefft_color(j, pid_m):
            if not match_colors:
                return None
            if match_by == "id" and pid_m is not None:
                return _color_for_id(("pair", pid_m))
            if match_by == "position" and j < len(tefft_lines):
                return tefft_lines[j][0].get_color()
            return None

        # -------- 
        # Plot MUELLER (dotted) 
        # --------
        for j, p in enumerate(m_list):
            if p is None:
                continue
            age = p["age_years"]
            # Allow vehicle_type -> HLE mapping if hle_cm not provided
            if "hle_cm" in p and p["hle_cm"] is not None:
                hle = p["hle_cm"]
            else:
                vtype = p.get("vehicle_type", None)
                if vtype is None:
                    raise ValueError("Mueller profile needs 'hle_cm' or 'vehicle_type'.")
                hle = vehicle_hle_guess(vtype, suv_as=p.get("suv_as", "pickup"), van_as=p.get("van_as", "pickup"))
            sm  = p.get("sex_male", 0)
            sev = p.get("severity", mueller_severity)
            pid = p.get(id_key, None)

            y = mueller_injury_probability(x_kmh, hle, sm, age, severity=sev)
            color = _paired_tefft_color(j, pid)

            (line,) = ax.plot(x, y, linewidth=2.6, linestyle=':',
                              label=_label_mueller(age, hle, sm, sev, pid, include_id_in_label),
                              color=color)

            for q in q_levels:
                xq = _x_at_y(x, y, q)
                if xq is None:
                    continue
                ax.scatter([xq], [q], s=marker_size, color=line.get_color(), **marker_kwargs)
                if show_marker_labels:
                    angle = _resolve_rotation(marker_label_rotation, x, y, xq)
                    ann = ax.annotate(
                        f"  {int(round(q*100))}% @ {xq:.1f} {units}",
                        xy=(xq, q), xytext=marker_label_offset, textcoords="offset points",
                        va="bottom", ha="left", fontsize=9,
                        rotation=angle, rotation_mode="anchor",
                    )
                    if use_text_stroke:
                        ann.set_path_effects([pe.withStroke(linewidth=3, foreground="black", alpha=0.6)])

        # -------- 
        # Axes & output 
        # --------
        ax.set_ylim(0, 1)
        ax.set_xlabel(x_label)
        ax.set_ylabel("Probability (death/injury)")
        ax.grid(True, alpha=0.35)
        if title is None:
            if model == "both":
                the_title = "Pedestrian risk vs. speed — Tefft (2013) & Monfort/Mueller (2025)"
            elif model == "tefft":
                the_title = "Pedestrian death risk vs. speed — Tefft (2013)"
            else:
                the_title = f"Pedestrian injury risk vs. speed — Mueller (2025) [{mueller_severity}]"
        else:
            the_title = title
        ax.set_title(the_title)
        ax.legend()

        if out:
            plt.savefig(out, bbox_inches="tight", dpi=160)
            print(f"Saved figure to: {out}")
        if show and created_fig:
            plt.show()

    return x, ax








