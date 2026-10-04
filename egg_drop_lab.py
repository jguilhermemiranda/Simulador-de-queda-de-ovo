"""
Egg Drop Lab 3.4
================
Simulação de queda, dimensionamento de paraquedas e proteção para um ovo.

Mudanças em relação à 3.3 (resumo; detalhes no relatório de bugs):
  - collect_data() quebrava com KeyError (vento) -> corrigido.
  - Clique nos gráficos quebrava (event.inaxes em evento Tk) -> corrigido.
  - Validação numérica estrita (antes, texto inválido virava valor padrão).
  - Resultado de simulação "velho" não é mais reutilizado pela Proteção.
  - Interpolação do instante de impacto no integrador RK4.
  - Interface redesenhada na paleta roxa do projeto, com abas roláveis.
  - Marca d'água / carimbo I.D.Y GUY (módulo brand.py).
"""
import csv
import json
import math
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from brand import BRAND


APP_TITLE = "Egg Drop Lab 3.4"

# --------------------------------------------------------------------------
# Paleta (preferência do projeto)
# --------------------------------------------------------------------------
COLOR1 = "#6a2f68"
COLOR2 = "#50234e"
COLOR3 = "#351834"
COLOR4 = "#1b0c1a"
COLOR5 = "#000000"

BG = COLOR5
PANEL = COLOR4
PANEL_2 = COLOR3
BORDER = COLOR2
ACCENT = COLOR1

# Tons derivados da paleta, necessários só para legibilidade de texto sobre
# fundos escuros (COLOR1 sobre COLOR4 tem contraste baixo demais para texto).
ACCENT_HOVER = "#85408a"
ACCENT_LIGHT = "#d9a8d6"
DIM_COLOR = "#f0c8ee"
TEXT = "#f3e8f2"
MUTED = "#b79fb5"

# Cores semânticas (status). Mantidas fora da paleta de propósito: verde /
# amarelo / vermelho precisam ser reconhecíveis como aprovado / atenção / falha.
GREEN = "#34d399"
YELLOW = "#fbbf24"
RED = "#fb7185"

FONT = ("Segoe UI", 10)
MONO = ("Consolas", 10)

CROSS_AREA_FACTOR = 0.60   # área do paraquedas cruciforme ≈ 0,60 · env · alt
MAX_SIM_TIME = 120.0

# Parâmetros que ficam travados no modo Básico.
ADVANCED_KEYS = (
    "gravity",
    "rho",
    "cd",
    "safety_factor",
    "para_cd",
    "structural_stop_mm",
    "crush_stress_kpa",
)

# Valores padrão: fonte única (antes estavam duplicados em add_field e
# _set_defaults).
DEFAULTS = {
    "height_m": 3.0,
    "gravity": 9.80665,
    "rho": 1.225,
    "wind_speed_ms": 0.0,
    "wind_uncertainty_pct": 20,
    "egg_mass_g": 60,
    "egg_length_mm": 60,
    "egg_width_mm": 45,
    "egg_height_mm": 45,
    "frontal_area_cm2": 20,
    "shell_strength_kpa": 300,
    "allowed_g": 100,
    "structure_mass_g": 100,
    "container_length_cm": 15,
    "container_width_cm": 15,
    "container_height_cm": 15,
    "cd": 0.47,
    "safety_factor": 2.0,
    "para_enabled": "Não",
    "parachute_type": "Redondo",
    "target_velocity_ms": 4.0,
    "deploy_height_m": 1.5,
    "para_cd": 1.5,
    "max_diameter_cm": 60,
    "parachute_width_cm": 50,
    "parachute_height_cm": 50,
    "cord_length_cm": 45,
    "cord_count": "4",
    "padding_mm": 40,
    "compression": 70,
    "structural_stop_mm": 5,
    "crush_stress_kpa": 80,
    "protection_area_cm2": 30,
    "target_margin": 1.20,
    "opt_min_padding_mm": 5,
    "opt_max_padding_mm": 80,
    "opt_step_mm": 5,
    "opt_max_structure_g": 300,
}


# ==========================================================================
# Utilidades
# ==========================================================================
def clamp(x, lo, hi):
    return max(lo, min(hi, x))


def safe_float(value, default=None):
    try:
        return float(str(value).replace(",", "."))
    except (ValueError, TypeError):
        return default


def fmt(x, digits=2):
    if x is None or not math.isfinite(x):
        return "-"
    return f"{x:.{digits}f}".replace(".", ",")


def fmt_si(x):
    if x is None or not math.isfinite(x):
        return "-"
    ax = abs(x)
    if ax >= 1000:
        return f"{x/1000:.2f} k".replace(".", ",")
    if ax < 0.001 and ax != 0:
        return f"{x*1000:.2f} m".replace(".", ",")
    return f"{x:.3f}".replace(".", ",")


# ==========================================================================
# Física
# ==========================================================================
def simulate_fall(
    height,
    mass,
    area,
    cd,
    rho,
    gravity,
    parachute_area=0.0,
    parachute_cd=0.0,
    parachute_deploy_height=0.0,
    dt=0.002,
    max_time=MAX_SIM_TIME,
):
    """
    Simulação RK4. y = altura acima do solo, v = velocidade para baixo.
    O instante e a velocidade de impacto são interpolados dentro do último
    passo (antes o passo inteiro era executado e y era truncado em 0,
    superestimando tempo e velocidade de impacto).
    """
    if min(height, mass, area, gravity, cd, rho) <= 0:
        raise ValueError(
            "Altura, massa, área, gravidade, Cd e densidade devem ser positivos."
        )

    has_para = parachute_area > 0 and parachute_cd > 0

    def acceleration(y_, v_):
        cda = cd * area
        if has_para and y_ <= parachute_deploy_height:
            cda += parachute_cd * parachute_area
        drag = 0.5 * rho * cda * abs(v_) * v_
        return gravity - drag / mass

    t, y, v = 0.0, height, 0.0
    ts, ys, vs = [t], [y], [v]
    accs, energies = [gravity], [0.0]
    opening_time = None
    landed = False

    while t < max_time:
        if has_para and opening_time is None and y <= parachute_deploy_height:
            opening_time = t

        k1y, k1v = -v, acceleration(y, v)
        k2y = -(v + 0.5 * dt * k1v)
        k2v = acceleration(y + 0.5 * dt * k1y, v + 0.5 * dt * k1v)
        k3y = -(v + 0.5 * dt * k2v)
        k3v = acceleration(y + 0.5 * dt * k2y, v + 0.5 * dt * k2v)
        k4y = -(v + dt * k3v)
        k4v = acceleration(y + dt * k3y, v + dt * k3v)

        y_new = y + dt * (k1y + 2 * k2y + 2 * k3y + k4y) / 6
        v_new = v + dt * (k1v + 2 * k2v + 2 * k3v + k4v) / 6

        if y_new <= 0:
            frac = y / (y - y_new)          # y > 0 e y_new <= 0 -> denominador > 0
            t_hit = t + frac * dt
            v_hit = v + frac * (v_new - v)
            ts.append(t_hit)
            ys.append(0.0)
            vs.append(v_hit)
            accs.append(acceleration(0.0, v_hit))
            energies.append(0.5 * mass * v_hit * v_hit)
            landed = True
            break

        y, v, t = y_new, v_new, t + dt
        ts.append(t)
        ys.append(y)
        vs.append(v)
        accs.append(acceleration(y, v))
        energies.append(0.5 * mass * v * v)

    vt_free = math.sqrt(2 * mass * gravity / (rho * cd * area))
    vt_final = None
    if has_para:
        vt_final = math.sqrt(
            2 * mass * gravity
            / (rho * (cd * area + parachute_cd * parachute_area))
        )

    return {
        "time": np.array(ts),
        "height": np.array(ys),
        "velocity": np.array(vs),
        "acceleration": np.array(accs),
        "energy": np.array(energies),
        "fall_time": ts[-1],
        "impact_velocity": vs[-1],
        "opening_time": opening_time,
        "landed": landed,
        "terminal_velocity": vt_free,          # só corpo (igual à 3.3)
        "terminal_velocity_final": vt_final,   # corpo + paraquedas (novo)
    }


def parachute_area_for_velocity(mass, gravity, rho, cd, target_velocity):
    if min(mass, gravity, rho, cd, target_velocity) <= 0:
        raise ValueError("Todos os parâmetros do paraquedas devem ser positivos.")
    return 2 * mass * gravity / (rho * cd * target_velocity**2)


def circular_diameter(area):
    return 2 * math.sqrt(area / math.pi)


def parachute_dimensions(ptype, area, ratio=1.0):
    """Retorna (largura, altura) em metros para a área dada."""
    if ptype == "Redondo":
        d = circular_diameter(area)
        return d, d
    ratio = max(ratio, 0.1)
    factor = 1.0 if ptype == "Retangular" else CROSS_AREA_FACTOR
    h = math.sqrt(area / (ratio * factor))
    return ratio * h, h


def parachute_area_from_size(ptype, max_dim, ratio=1.0):
    """Área para um paraquedas cuja MAIOR dimensão é max_dim (m)."""
    if ptype == "Redondo":
        return math.pi * (max_dim / 2) ** 2
    ratio = max(ratio, 0.1)
    if ratio >= 1:
        w, h = max_dim, max_dim / ratio
    else:
        w, h = max_dim * ratio, max_dim
    factor = 1.0 if ptype == "Retangular" else CROSS_AREA_FACTOR
    return factor * w * h


def impact_analysis(
    mass,
    impact_velocity,
    stop_distance,
    frontal_area,
    shell_strength,
    safety_factor,
    allowed_g,
    gravity,
):
    stop_distance = max(stop_distance, 1e-5)
    frontal_area = max(frontal_area, 1e-8)

    energy = 0.5 * mass * impact_velocity**2

    # Desaceleração cinemática líquida (aceleração constante na parada).
    net_decel = impact_velocity**2 / (2 * stop_distance)

    # Reação estrutural: o peso continua atuando durante a parada.
    decel = net_decel + gravity
    g_load = decel / gravity
    force = mass * decel
    pressure = force / frontal_area

    allowable_force = shell_strength * frontal_area / max(safety_factor, 1e-9)
    force_margin = allowable_force / max(force, 1e-9)
    g_margin = allowed_g / max(g_load, 1e-9)

    force_score = force / max(allowable_force, 1e-9)
    g_score = g_load / max(allowed_g, 1e-9)
    heuristic_index = clamp(100 * (0.6 * force_score + 0.4 * g_score), 0, 100)

    if heuristic_index < 20:
        label = "MUITO BAIXO"
    elif heuristic_index < 40:
        label = "BAIXO"
    elif heuristic_index < 60:
        label = "MODERADO"
    elif heuristic_index < 80:
        label = "ALTO"
    else:
        label = "MUITO ALTO"

    return {
        "impact_energy": energy,
        "stop_distance": stop_distance,
        "deceleration": decel,
        "net_deceleration": net_decel,
        "g_load": g_load,
        "force": force,
        "pressure": pressure,
        "allowable_force": allowable_force,
        "force_margin": force_margin,
        "g_margin": g_margin,
        "heuristic_index": heuristic_index,
        "risk_label": label,
    }


def energy_absorption_capacity(padding, compression, crush_stress, area):
    """E_abs ≈ σ_crush · A · d_compressão (estimativa simplificada)."""
    compression_distance = max(0.0, padding * clamp(compression, 0, 1))
    force = max(0.0, crush_stress) * max(0.0, area)
    return force * compression_distance, force, compression_distance


def estimate_landing_position(d, sim):
    fall_time = max(0.0, float(sim["fall_time"]))
    wind = max(0.0, float(d.get("wind_speed", 0.0)))
    uncertainty = max(0.0, float(d.get("wind_uncertainty", 0.20)))

    # Modelo simplificado: deriva = vento · tempo de queda · acoplamento.
    # Fornece uma faixa de projeto, não uma previsão meteorológica.
    coupling = 1.0 if d.get("para_enabled") else 0.55
    position = wind * fall_time * coupling
    spread = max(0.02, position * uncertainty)
    return {
        "position_m": position,
        "min_m": max(0.0, position - spread),
        "max_m": position + spread,
    }


def compute_result(d):
    """Simulação + impacto + proteção, sem tocar na interface."""
    parachute_area = 0.0
    if d["para_enabled"]:
        parachute_area = parachute_area_for_velocity(
            d["total_mass"], d["gravity"], d["rho"],
            d["para_cd"], d["target_velocity"],
        )

    sim = simulate_fall(
        height=d["height"],
        mass=d["total_mass"],
        area=d["frontal_area"],
        cd=d["cd"],
        rho=d["rho"],
        gravity=d["gravity"],
        parachute_area=parachute_area,
        parachute_cd=d["para_cd"],
        parachute_deploy_height=d["deploy_height"],
    )

    stop_distance = d["padding"] * d["compression"] + d["structural_stop"]

    impact = impact_analysis(
        mass=d["total_mass"],
        impact_velocity=sim["impact_velocity"],
        stop_distance=stop_distance,
        frontal_area=d["frontal_area"],
        shell_strength=d["shell_strength"],
        safety_factor=d["safety_factor"],
        allowed_g=d["allowed_g"],
        gravity=d["gravity"],
    )

    capacity, _, _ = energy_absorption_capacity(
        d["padding"], d["compression"], d["crush_stress"], d["protection_area"]
    )

    geometry_ok = (
        d["egg_length"] + 2 * d["padding"] <= d["container_length"]
        and d["egg_width"] + 2 * d["padding"] <= d["container_width"]
        and d["egg_height"] + 2 * d["padding"] <= d["container_height"]
    )

    return {
        "data": d,
        "landing": estimate_landing_position(d, sim),
        "simulation": sim,
        "impact": impact,
        "absorption_capacity": capacity,
        "energy_margin": capacity / max(impact["impact_energy"], 1e-9),
        "geometry_ok": geometry_ok,
    }


# ==========================================================================
# Desenho (funções puras: recebem uma Figure)
# ==========================================================================
def style_axis(ax):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=MUTED)
    for spine in ax.spines.values():
        spine.set_color(BORDER)
    ax.grid(alpha=0.35, color=BORDER)


GRAPH_TITLES = [
    "Velocidade",
    "Aceleração",
    "Altura",
    "Energia cinética",
    "Velocidade × altura",
    "Previsão de aterrissagem horizontal",
]
_SERIES = [
    ("velocity", "m/s"),
    ("acceleration", "m/s²"),
    ("height", "m"),
    ("energy", "J"),
]


def draw_graph(ax, index, sim, landing, big=False):
    lw = 2.6 if big else 1.8
    style_axis(ax)

    if index < 4:
        key, unit = _SERIES[index]
        ax.plot(sim["time"], sim[key], color=ACCENT_LIGHT, linewidth=lw)
        ax.set_xlabel("Tempo (s)", color=MUTED)
        ax.set_ylabel(unit, color=MUTED)
        if sim["opening_time"] is not None:
            ax.axvline(sim["opening_time"], linestyle="--", linewidth=1.2, color=YELLOW)
            ax.text(sim["opening_time"], ax.get_ylim()[1], " abertura",
                    color=YELLOW, va="top")
    elif index == 4:
        ax.plot(sim["height"], sim["velocity"], color=ACCENT_LIGHT, linewidth=lw)
        ax.invert_xaxis()   # o ovo desce: altura diminui da esquerda p/ direita
        ax.set_xlabel("Altura (m)", color=MUTED)
        ax.set_ylabel("Velocidade (m/s)", color=MUTED)
    else:
        x, xmin, xmax = landing["position_m"], landing["min_m"], landing["max_m"]
        ax.axvspan(xmin, xmax, alpha=0.35, color=ACCENT)
        ax.axvline(x, linewidth=lw, color=ACCENT_LIGHT)
        ax.scatter([x], [0], s=110 if big else 60, zorder=5, color=ACCENT_LIGHT)
        ax.annotate(f"estimativa: {x:.2f} m", (x, 0), xytext=(x, 0.55),
                    color=TEXT, ha="center",
                    arrowprops=dict(arrowstyle="->", color=ACCENT_LIGHT))
        ax.annotate(f"faixa: {xmin:.2f}–{xmax:.2f} m", (x, 0), xytext=(x, 0.25),
                    color=MUTED, ha="center")
        ax.set_xlabel("Distância do ponto de lançamento (m)", color=MUTED)
        ax.set_xlim(0, max(1.0, xmax * 1.35))
        ax.set_ylim(-0.15, 0.8)
        ax.set_yticks([])

    ax.set_title(GRAPH_TITLES[index], color=TEXT, fontsize=16 if big else 11)


def draw_all_graphs(fig, sim, landing):
    fig.clear()
    fig.patch.set_facecolor(BG)
    axes = fig.subplots(2, 3)
    for i, ax in enumerate(axes.flat):
        draw_graph(ax, i, sim, landing)
    fig.tight_layout(pad=2)
    BRAND.watermark_figure(fig)


def draw_single_graph(fig, index, sim, landing):
    fig.clear()
    fig.patch.set_facecolor(BG)
    ax = fig.add_subplot(111)
    draw_graph(ax, index, sim, landing, big=True)
    fig.tight_layout(pad=2)
    BRAND.watermark_figure(fig)


def _dim_line(ax, p1, p2, text, label_offset=(0, 0), ha="center", va="bottom"):
    ax.annotate("", xy=p1, xytext=p2,
                arrowprops=dict(arrowstyle="<->", color=DIM_COLOR, lw=1.2))
    mx, my = (p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2
    ax.text(mx + label_offset[0], my + label_offset[1], text,
            color=TEXT, ha=ha, va=va, fontsize=9)


def _cross_arm_width(w, h):
    """Largura do braço tal que a área do desenho = CROSS_AREA_FACTOR·w·h."""
    s = w + h
    disc = max(s * s - 4 * CROSS_AREA_FACTOR * w * h, 0.0)
    return min((s - math.sqrt(disc)) / 2, min(w, h))


def _attachment_points(ptype, w, h, n):
    if ptype == "Redondo":
        return [((w / 2) * math.cos(a), (h / 2) * math.sin(a))
                for a in np.linspace(0, 2 * np.pi, n, endpoint=False)]
    if ptype == "Retangular":
        pts = [(-w/2, -h/2), (w/2, -h/2), (w/2, h/2), (-w/2, h/2),
               (0, -h/2), (0, h/2), (-w/2, 0), (w/2, 0)]
        return pts[:n]
    a = _cross_arm_width(w, h)
    pts = [(w/2, 0), (0, h/2), (-w/2, 0), (0, -h/2),
           (a/2, a/2), (-a/2, a/2), (-a/2, -a/2), (a/2, -a/2)]
    return pts[:n]


def draw_parachute(fig, ptype, w, h, cord_count, cord_length):
    fig.clear()
    fig.patch.set_facecolor(BG)
    ax_top, ax_side = fig.subplots(1, 2)
    style_axis(ax_top)
    style_axis(ax_side)

    # ---- vista superior ------------------------------------------------
    if ptype == "Redondo":
        th = np.linspace(0, 2 * np.pi, 300)
        xs, ys = (w / 2) * np.cos(th), (h / 2) * np.sin(th)
    elif ptype == "Retangular":
        xs = np.array([-w/2, w/2, w/2, -w/2, -w/2])
        ys = np.array([-h/2, -h/2, h/2, h/2, -h/2])
    else:
        a = _cross_arm_width(w, h)
        poly = [(-a/2, -h/2), (a/2, -h/2), (a/2, -a/2), (w/2, -a/2),
                (w/2, a/2), (a/2, a/2), (a/2, h/2), (-a/2, h/2),
                (-a/2, a/2), (-w/2, a/2), (-w/2, -a/2), (-a/2, -a/2),
                (-a/2, -h/2)]
        xs, ys = np.array([p[0] for p in poly]), np.array([p[1] for p in poly])

    ax_top.fill(xs, ys, color=ACCENT, alpha=0.35)
    ax_top.plot(xs, ys, color=ACCENT_LIGHT, linewidth=2.2)

    for px, py in _attachment_points(ptype, w, h, cord_count):
        ax_top.plot([px, 0], [py, 0], color=YELLOW, alpha=0.5, linewidth=1)
        ax_top.scatter([px], [py], s=40, color=YELLOW, zorder=5)
    ax_top.scatter([0], [0], s=60, color=DIM_COLOR, zorder=6)
    ax_top.text(0, 0, " carga", color=TEXT, fontsize=9, va="bottom")

    m = max(w, h) * 0.18
    if ptype == "Redondo":
        _dim_line(ax_top, (-w/2, -h/2 - m), (w/2, -h/2 - m),
                  f"Ø {w*100:.1f} cm", (0, -m * 0.2), va="top")
    else:
        _dim_line(ax_top, (-w/2, -h/2 - m), (w/2, -h/2 - m),
                  f"{w*100:.1f} cm", (0, -m * 0.2), va="top")
        _dim_line(ax_top, (w/2 + m, -h/2), (w/2 + m, h/2),
                  f"{h*100:.1f} cm", (m * 0.2, 0), ha="left", va="center")
    pad = max(w, h) * 0.45
    ax_top.set_xlim(-w/2 - pad, w/2 + pad)
    ax_top.set_ylim(-h/2 - pad, h/2 + pad)
    ax_top.set_aspect("equal", adjustable="box")
    ax_top.set_title("VISTA SUPERIOR — cotas e fixação", color=TEXT)
    ax_top.set_xlabel("pontos amarelos = fixação das cordas", color=MUTED)

    # ---- vista lateral (em escala) -----------------------------------------
    ax_side.plot([-w/2, w/2], [0, 0], color=ACCENT_LIGHT, linewidth=3)
    for px in np.linspace(-w/2, w/2, max(2, cord_count)):
        ax_side.plot([px, 0], [0, -cord_length], color=YELLOW, alpha=0.65, linewidth=1.2)
        ax_side.scatter([px], [0], s=36, color=YELLOW, zorder=5)
    ax_side.scatter([0], [-cord_length], s=80, color=DIM_COLOR, zorder=5)
    ax_side.text(0.02 * w, -cord_length, " carga", color=TEXT, va="center")

    big = max(w, cord_length)
    _dim_line(ax_side, (w/2 + big * 0.12, 0), (w/2 + big * 0.12, -cord_length),
              f"{cord_length*100:.1f} cm", (big * 0.02, 0), ha="left", va="center")
    _dim_line(ax_side, (-w/2, big * 0.10), (w/2, big * 0.10),
              f"{w*100:.1f} cm", (0, big * 0.01))
    ax_side.set_xlim(-big * 0.75, big * 0.75)
    ax_side.set_ylim(-cord_length - big * 0.15, big * 0.30)
    ax_side.set_aspect("equal", adjustable="box")
    ax_side.set_title("VISTA LATERAL — comprimento das cordas", color=TEXT)
    ax_side.set_xlabel("desenho em escala", color=MUTED)

    fig.tight_layout(pad=2)
    BRAND.watermark_figure(fig)


# ==========================================================================
# Widgets auxiliares
# ==========================================================================
class ScrollFrame(ttk.Frame):
    """Frame com rolagem vertical; coloque o conteúdo em .body."""

    def __init__(self, parent):
        super().__init__(parent)
        self.canvas = tk.Canvas(self, bg=BG, highlightthickness=0)
        self.vsb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.body = ttk.Frame(self.canvas)
        self._win = self.canvas.create_window((0, 0), window=self.body, anchor="nw")
        self.canvas.configure(yscrollcommand=self.vsb.set)
        self.canvas.pack(side="left", fill="both", expand=True)
        self.vsb.pack(side="right", fill="y")
        self.body.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")),
        )
        self.canvas.bind(
            "<Configure>",
            lambda e: self.canvas.itemconfigure(self._win, width=e.width),
        )

    def scroll(self, units):
        bbox = self.canvas.bbox("all")
        if bbox and (bbox[3] - bbox[1]) > self.canvas.winfo_height():
            self.canvas.yview_scroll(units, "units")


def set_text(widget, text):
    """Escreve em um Text somente-leitura."""
    widget.configure(state="normal")
    widget.delete("1.0", "end")
    widget.insert("end", text)
    widget.configure(state="disabled")


def make_text(parent, height=10, font=MONO):
    t = tk.Text(
        parent, height=height, bg=PANEL_2, fg=TEXT, relief="flat",
        insertbackground=TEXT, selectbackground=ACCENT, font=font,
        padx=14, pady=12, wrap="word", state="disabled",
    )
    return t


def configure_style(root):
    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except tk.TclError:
        pass

    root.option_add("*TCombobox*Listbox.background", PANEL_2)
    root.option_add("*TCombobox*Listbox.foreground", TEXT)
    root.option_add("*TCombobox*Listbox.selectBackground", ACCENT)
    root.option_add("*TCombobox*Listbox.selectForeground", "#ffffff")

    style.configure(".", background=BG, foreground=TEXT, fieldbackground=PANEL_2,
                    font=FONT, bordercolor=BORDER, lightcolor=BORDER,
                    darkcolor=BORDER, troughcolor=BG)
    style.configure("TFrame", background=BG)
    style.configure("Panel.TFrame", background=PANEL)
    style.configure("TLabel", background=BG, foreground=TEXT)
    style.configure("Muted.TLabel", background=BG, foreground=MUTED)
    style.configure("Panel.TLabel", background=PANEL, foreground=TEXT)
    style.configure("PanelMuted.TLabel", background=PANEL, foreground=MUTED)
    style.configure("Title.TLabel", background=BG, foreground=TEXT,
                    font=("Segoe UI", 20, "bold"))
    style.configure("SubTitle.TLabel", background=BG, foreground=MUTED)
    style.configure("CardTitle.TLabel", background=PANEL, foreground=ACCENT_LIGHT,
                    font=("Segoe UI", 11, "bold"))
    style.configure("Footer.TLabel", background=PANEL, foreground=MUTED,
                    padding=(14, 6))

    style.configure("Accent.TButton", background=ACCENT, foreground="#ffffff",
                    borderwidth=0, padding=(14, 9), font=("Segoe UI", 10, "bold"),
                    lightcolor=ACCENT, darkcolor=ACCENT, bordercolor=ACCENT)
    style.map("Accent.TButton",
              background=[("active", ACCENT_HOVER), ("pressed", COLOR2)])
    style.configure("Secondary.TButton", background=PANEL_2, foreground=TEXT,
                    borderwidth=0, padding=(12, 8), lightcolor=PANEL_2,
                    darkcolor=PANEL_2, bordercolor=BORDER)
    style.map("Secondary.TButton", background=[("active", COLOR2)])

    style.configure("TEntry", fieldbackground=PANEL_2, foreground=TEXT,
                    insertcolor=TEXT, padding=6, bordercolor=BORDER,
                    lightcolor=PANEL_2, darkcolor=PANEL_2)
    style.map("TEntry",
              fieldbackground=[("disabled", PANEL)],
              foreground=[("disabled", MUTED)],
              bordercolor=[("focus", ACCENT_LIGHT)])
    style.configure("TCombobox", fieldbackground=PANEL_2, foreground=TEXT,
                    background=PANEL_2, arrowcolor=TEXT, padding=5,
                    selectbackground=PANEL_2, selectforeground=TEXT,
                    bordercolor=BORDER, lightcolor=PANEL_2, darkcolor=PANEL_2)
    style.map("TCombobox",
              fieldbackground=[("disabled", PANEL), ("readonly", PANEL_2)],
              foreground=[("disabled", MUTED)],
              background=[("active", COLOR2)])

    style.configure("TNotebook", background=BG, borderwidth=0)
    style.configure("TNotebook.Tab", background=PANEL, foreground=MUTED,
                    padding=(14, 8), borderwidth=0)
    style.map("TNotebook.Tab",
              background=[("selected", ACCENT), ("active", COLOR2)],
              foreground=[("selected", "#ffffff")])

    style.configure("Treeview", background=PANEL_2, fieldbackground=PANEL_2,
                    foreground=TEXT, rowheight=26, borderwidth=0)
    style.configure("Treeview.Heading", background=COLOR2, foreground=TEXT,
                    relief="flat", padding=6, font=("Segoe UI", 10, "bold"))
    style.map("Treeview", background=[("selected", ACCENT)],
              foreground=[("selected", "#ffffff")])
    style.map("Treeview.Heading", background=[("active", ACCENT)])

    style.configure("Vertical.TScrollbar", background=COLOR2, troughcolor=BG,
                    arrowcolor=TEXT, bordercolor=BG)


# ==========================================================================
# Aplicação
# ==========================================================================
class EggDropLab:
    def __init__(self, root):
        self.root = root
        root.title(APP_TITLE)
        root.geometry("1420x900")
        root.minsize(1180, 760)
        try:
            root.state("zoomed")          # Windows: maximizado, com barra de título
        except tk.TclError:
            pass
        root.bind("<F11>", self._toggle_fullscreen)
        root.bind("<Escape>", self._exit_fullscreen)
        root.configure(bg=BG)

        self.vars = {}
        self.widgets = {}
        self.labels = {}
        self.label_text = {}
        self.tabs = {}
        self.kpis = {}
        self.result = None
        self.last_project_path = None
        self.mode = tk.StringVar(value="Básico")
        self._last_mode = "Básico"
        self.status_var = tk.StringVar(value="Pronto.")

        configure_style(root)
        BRAND.set_window_icon(root)
        self._build_header()
        self._build_footer()
        self._build_main()
        self._set_defaults()
        self.refresh_parachute_drawing()

        root.bind_all("<MouseWheel>", self._on_mousewheel)
        root.bind_all("<Button-4>", self._on_mousewheel)
        root.bind_all("<Button-5>", self._on_mousewheel)

    # ------------------------------------------------------------------
    # Janela
    # ------------------------------------------------------------------
    def _toggle_fullscreen(self, event=None):
        self.root.attributes("-fullscreen", not bool(self.root.attributes("-fullscreen")))

    def _exit_fullscreen(self, event=None):
        self.root.attributes("-fullscreen", False)

    def _on_mousewheel(self, event):
        if getattr(event, "num", None) == 4:
            units = -2
        elif getattr(event, "num", None) == 5:
            units = 2
        else:
            units = -2 if event.delta > 0 else 2
        try:
            w = self.root.winfo_containing(event.x_root, event.y_root)
        except (KeyError, tk.TclError):
            return
        while w is not None:
            if isinstance(w, tk.Text):
                return
            if isinstance(w, ScrollFrame):
                w.scroll(units)
                return
            w = getattr(w, "master", None)

    def set_status(self, text, color=MUTED):
        self.status_var.set(text)
        self.footer.configure(foreground=color)

    # ------------------------------------------------------------------
    # Construção base
    # ------------------------------------------------------------------
    def _build_header(self):
        header = ttk.Frame(self.root)
        header.pack(fill="x", padx=20, pady=(16, 8))

        logo = BRAND.tk_image(56, background=BG)
        if logo:
            ttk.Label(header, image=logo).pack(side="left", padx=(0, 12))

        left = ttk.Frame(header)
        left.pack(side="left")
        ttk.Label(left, text="🥚 EGG DROP LAB", style="Title.TLabel").pack(anchor="w")
        ttk.Label(left, text="Simulação, paraquedas e dimensionamento de proteção",
                  style="SubTitle.TLabel").pack(anchor="w")

        right = ttk.Frame(header)
        right.pack(side="right")
        ttk.Label(right, text="Modo:", style="Muted.TLabel").pack(side="left", padx=(0, 6))
        combo = ttk.Combobox(right, textvariable=self.mode,
                             values=["Básico", "Avançado"], state="readonly", width=12)
        combo.pack(side="left")
        combo.bind("<<ComboboxSelected>>", lambda _: self._on_mode_change())

        tk.Frame(self.root, bg=ACCENT, height=2).pack(fill="x", padx=20)

    def _build_footer(self):
        self.footer = ttk.Label(self.root, textvariable=self.status_var,
                                style="Footer.TLabel", anchor="w")
        self.footer.pack(side="bottom", fill="x")

    def _build_main(self):
        container = ttk.Frame(self.root)
        container.pack(fill="both", expand=True, padx=20, pady=(10, 10))

        self.notebook = ttk.Notebook(container)
        self.notebook.pack(fill="both", expand=True)

        self.tab_home = self._add_tab("home", "🏠 Início")
        self.tab_project = self._add_tab("project", "🥚 Parâmetros")
        self.tab_para = self._add_tab("para", "☂️ Paraquedas")
        self.tab_protection = self._add_tab("protection", "🧱 Proteção")
        self.tab_sim = self._add_tab("sim", "🧪 Simulação")
        self.tab_graphs = self._add_tab("graphs", "📊 Gráficos", scroll=False)
        self.tab_equations = self._add_tab("equations", "🧮 Equações", scroll=False)
        self.tab_save = self._add_tab("save", "💾 Arquivo e otimização")

        self._build_home()
        self._build_project()
        self._build_para()
        self._build_protection()
        self._build_sim()
        self._build_graphs()
        self._build_equations()
        self._build_save()

    def _add_tab(self, name, text, scroll=True):
        outer = ttk.Frame(self.notebook)
        self.notebook.add(outer, text=text)
        self.tabs[name] = outer
        if not scroll:
            return outer
        sf = ScrollFrame(outer)
        sf.pack(fill="both", expand=True)
        return sf.body

    def goto(self, name):
        self.notebook.select(self.tabs[name])

    def card(self, parent, title, row, column=0, columnspan=1):
        outer = tk.Frame(parent, bg=PANEL, highlightbackground=BORDER,
                         highlightthickness=1)
        outer.grid(row=row, column=column, columnspan=columnspan,
                   sticky="nsew", padx=7, pady=7)
        inner = ttk.Frame(outer, style="Panel.TFrame", padding=16)
        inner.pack(fill="both", expand=True)
        ttk.Label(inner, text=title, style="CardTitle.TLabel").pack(anchor="w", pady=(0, 10))
        return inner

    def add_field(self, parent, key, label, unit="", values=None):
        line = ttk.Frame(parent, style="Panel.TFrame")
        line.pack(fill="x", pady=4)

        lab = ttk.Label(line, text=label, style="Panel.TLabel")
        lab.pack(side="left", fill="x", expand=True)

        var = tk.StringVar(value=str(DEFAULTS[key]))
        self.vars[key] = var

        if values:
            widget = ttk.Combobox(line, textvariable=var, values=values,
                                  state="readonly", width=13)
        else:
            widget = ttk.Entry(line, textvariable=var, width=14)
        widget.pack(side="right")

        if unit:
            ttk.Label(line, text=unit, style="PanelMuted.TLabel",
                      width=8, anchor="w").pack(side="right", padx=(7, 0))

        self.widgets[key] = widget
        self.labels[key] = lab
        self.label_text[key] = label
        return widget

    def _result_text(self, parent, height):
        t = make_text(parent, height=height)
        t.pack(fill="both", expand=True)
        return t

    # ------------------------------------------------------------------
    # Abas
    # ------------------------------------------------------------------
    def _build_home(self):
        t = self.tab_home
        t.columnconfigure(0, weight=1)

        hero = tk.Frame(t, bg=PANEL, highlightbackground=BORDER, highlightthickness=1)
        hero.grid(row=0, column=0, sticky="ew", padx=7, pady=7)
        inner = ttk.Frame(hero, style="Panel.TFrame", padding=28)
        inner.pack(fill="both", expand=True)

        mark = BRAND.tk_image(100, color=ACCENT, opacity=0.8, background=PANEL)
        if mark:
            tk.Label(inner, image=mark, bg=PANEL, bd=0).place(
                relx=1.0, x=-4, y=4, anchor="ne")

        ttk.Label(inner, text="🥚 Egg Drop Lab", style="Panel.TLabel",
                  font=("Segoe UI", 28, "bold")).pack(anchor="w")
        ttk.Label(inner, style="PanelMuted.TLabel", font=("Segoe UI", 11),
                  text="Projete, simule e compare soluções para proteger um ovo durante uma queda.",
                  ).pack(anchor="w", pady=(4, 15))

        buttons = ttk.Frame(inner, style="Panel.TFrame")
        buttons.pack(fill="x")
        actions = [
            ("🚀 Simular projeto", lambda: self.goto("sim")),
            ("☂️ Calcular paraquedas", lambda: self.goto("para")),
            ("🧱 Dimensionar proteção", lambda: self.goto("protection")),
            ("📊 Ver gráficos", lambda: self.goto("graphs")),
        ]
        for i, (text, cmd) in enumerate(actions):
            ttk.Button(buttons, text=text, command=cmd, style="Accent.TButton").grid(
                row=0, column=i, sticky="ew", padx=4)
            buttons.columnconfigure(i, weight=1)

        info = self.card(t, "Como funciona", 1)
        ttk.Label(
            info, style="Panel.TLabel", justify="left",
            text=(
                "1. Informe os dados do ovo e da estrutura na aba Parâmetros.\n"
                "2. Se quiser, ative o paraquedas e peça o dimensionamento automático.\n"
                "3. Simule a queda.\n"
                "4. Veja velocidade, aceleração, energia e margem de segurança.\n\n"
                "⚠️ Os resultados são estimativas de engenharia. "
                "Eles não garantem que um ovo real sobreviverá ao teste."
            ),
        ).pack(anchor="w")

        self.home_status = ttk.Label(t, text="Nenhuma simulação executada.",
                                     style="Muted.TLabel", wraplength=1100)
        self.home_status.grid(row=2, column=0, sticky="w", padx=14, pady=10)

    def _build_project(self):
        t = self.tab_project
        t.columnconfigure(0, weight=1)
        t.columnconfigure(1, weight=1)

        env = self.card(t, "🌍 Queda e ambiente", 0, 0)
        self.add_field(env, "height_m", "Altura da queda", "m")
        self.add_field(env, "gravity", "Gravidade", "m/s²")
        self.add_field(env, "rho", "Densidade do ar", "kg/m³")
        self.add_field(env, "wind_speed_ms", "Vento horizontal", "m/s")
        self.add_field(env, "wind_uncertainty_pct", "Incerteza do vento", "%")

        egg = self.card(t, "🥚 Ovo", 0, 1)
        self.add_field(egg, "egg_mass_g", "Massa", "g")
        self.add_field(egg, "egg_length_mm", "Comprimento", "mm")
        self.add_field(egg, "egg_width_mm", "Largura", "mm")
        self.add_field(egg, "egg_height_mm", "Altura", "mm")
        self.add_field(egg, "frontal_area_cm2", "Área frontal", "cm²")
        self.add_field(egg, "shell_strength_kpa", "Resistência estimada", "kPa")
        self.add_field(egg, "allowed_g", "Aceleração limite", "g")

        st = self.card(t, "📦 Estrutura / recipiente", 1, 0)
        self.add_field(st, "structure_mass_g", "Massa da estrutura", "g")
        self.add_field(st, "container_length_cm", "Comprimento interno", "cm")
        self.add_field(st, "container_width_cm", "Largura interna", "cm")
        self.add_field(st, "container_height_cm", "Altura interna", "cm")

        aero = self.card(t, "🌬️ Aerodinâmica", 1, 1)
        self.add_field(aero, "cd", "Coeficiente de arrasto")
        self.add_field(aero, "safety_factor", "Fator de segurança")

        buttons = ttk.Frame(t)
        buttons.grid(row=2, column=0, columnspan=2, sticky="ew", pady=8)
        buttons.columnconfigure(0, weight=1)
        buttons.columnconfigure(1, weight=1)
        ttk.Button(buttons, text="🧮 Calcular / simular", style="Accent.TButton",
                   command=self.run_simulation).grid(row=0, column=0, sticky="ew", padx=7)
        ttk.Button(buttons, text="☂️ Dimensionar paraquedas", style="Secondary.TButton",
                   command=self.calculate_parachute).grid(row=0, column=1, sticky="ew", padx=7)

    def _build_para(self):
        t = self.tab_para
        t.columnconfigure(0, weight=1)
        t.columnconfigure(1, weight=2)

        p = self.card(t, "☂️ Calculadora de paraquedas", 0, 0)
        self.add_field(p, "para_enabled", "Usar paraquedas", "", ["Não", "Sim"])
        self.add_field(p, "parachute_type", "Tipo de paraquedas", "",
                       ["Redondo", "Retangular", "Cruciforme"])
        self.add_field(p, "target_velocity_ms", "Velocidade alvo", "m/s")
        self.add_field(p, "deploy_height_m", "Altura de abertura", "m")
        self.add_field(p, "para_cd", "Cd do paraquedas")
        self.add_field(p, "max_diameter_cm", "Dimensão máxima", "cm")
        self.add_field(p, "parachute_width_cm", "Proporção — largura", "(rel.)")
        self.add_field(p, "parachute_height_cm", "Proporção — altura", "(rel.)")
        self.add_field(p, "cord_length_cm", "Comprimento das cordas", "cm")
        self.add_field(p, "cord_count", "Número de cordas", "", ["3", "4", "6", "8"])
        self.widgets["parachute_type"].bind(
            "<<ComboboxSelected>>", lambda _: self._on_parachute_type_change())

        ttk.Button(p, text="☂️ CALCULAR E ATUALIZAR DESENHO", style="Accent.TButton",
                   command=self.calculate_parachute).pack(fill="x", pady=(10, 0))

        r = self.card(t, "📐 Resultado", 0, 1)
        self.para_result = self._result_text(r, 14)

        dg = self.card(t, "📏 Projeto físico — cotas, vistas e pontos de fixação "
                          "(clique para ampliar)", 1, 0, 2)
        self.para_figure = Figure(figsize=(12, 4.6), dpi=100)
        self.para_canvas = FigureCanvasTkAgg(self.para_figure, master=dg)
        self.para_canvas.get_tk_widget().pack(fill="both", expand=True)
        self.para_canvas.get_tk_widget().bind(
            "<Button-1>", lambda e: self._open_parachute_fullscreen())

        tf = self.card(t, "🔎 Alternativas de dimensionamento", 2, 0, 2)
        cols = ("dimension", "area", "velocity", "energy")
        heads = {"dimension": "Maior dimensão", "area": "Área",
                 "velocity": "Velocidade terminal (só paraquedas)",
                 "energy": "Energia de impacto"}
        self.para_tree = ttk.Treeview(tf, columns=cols, show="headings", height=7)
        for c in cols:
            self.para_tree.heading(c, text=heads[c])
            self.para_tree.column(c, width=180, anchor="center")
        self.para_tree.tag_configure("ok", foreground=GREEN)
        self.para_tree.pack(fill="both", expand=True)

    def _build_protection(self):
        t = self.tab_protection
        t.columnconfigure(0, weight=1)
        t.columnconfigure(1, weight=1)

        prot = self.card(t, "🧱 Amortecimento", 0, 0)
        self.add_field(prot, "padding_mm", "Espessura da proteção", "mm")
        self.add_field(prot, "compression", "Compressão utilizável", "%")
        self.add_field(prot, "structural_stop_mm", "Parada estrutural", "mm")
        self.add_field(prot, "crush_stress_kpa", "Tensão de esmagamento", "kPa")
        self.add_field(prot, "protection_area_cm2", "Área efetiva", "cm²")

        lim = self.card(t, "🎯 Limites e resultado", 0, 1)
        self.add_field(lim, "target_margin", "Margem energética desejada", "x")
        self.protection_result = self._result_text(lim, 11)

        ttk.Button(t, text="🧱 DIMENSIONAR PROTEÇÃO", style="Accent.TButton",
                   command=self.calculate_protection).grid(
            row=1, column=0, columnspan=2, sticky="ew", padx=7, pady=8)

    def _build_sim(self):
        t = self.tab_sim
        t.columnconfigure(0, weight=1)

        self.sim_status_var = tk.StringVar(value="Nenhuma simulação executada.")
        self.sim_status_label = tk.Label(
            t, textvariable=self.sim_status_var, bg=PANEL, fg=MUTED,
            font=("Segoe UI", 12, "bold"), anchor="w", padx=16, pady=10,
            highlightbackground=BORDER, highlightthickness=1)
        self.sim_status_label.grid(row=0, column=0, sticky="ew", padx=7, pady=(7, 0))

        kpi_row = ttk.Frame(t)
        kpi_row.grid(row=1, column=0, sticky="ew")
        for i, (key, title) in enumerate([
            ("v", "Velocidade de impacto"), ("g", "Carga em g"),
            ("margin", "Margem energética"), ("risk", "Índice heurístico"),
        ]):
            f = tk.Frame(kpi_row, bg=PANEL, highlightbackground=BORDER,
                         highlightthickness=1)
            f.grid(row=0, column=i, sticky="nsew", padx=7, pady=7)
            kpi_row.columnconfigure(i, weight=1)
            tk.Label(f, text=title, bg=PANEL, fg=MUTED,
                     font=("Segoe UI", 9)).pack(anchor="w", padx=14, pady=(10, 0))
            val = tk.Label(f, text="—", bg=PANEL, fg=TEXT,
                           font=("Segoe UI", 20, "bold"))
            val.pack(anchor="w", padx=14, pady=(0, 10))
            self.kpis[key] = val

        summary = self.card(t, "📋 Resultado da simulação", 2)
        self.result_text = self._result_text(summary, 20)

        actions = ttk.Frame(t)
        actions.grid(row=3, column=0, sticky="ew", pady=7)
        for i, (text, cmd, style) in enumerate([
            ("🧪 SIMULAR", self.run_simulation, "Accent.TButton"),
            ("🔎 OTIMIZAR", self.optimize, "Secondary.TButton"),
            ("📊 VER GRÁFICOS", lambda: self.goto("graphs"), "Secondary.TButton"),
        ]):
            ttk.Button(actions, text=text, command=cmd, style=style).grid(
                row=0, column=i, sticky="ew", padx=7)
            actions.columnconfigure(i, weight=1)

    def _build_graphs(self):
        frame = ttk.Frame(self.tab_graphs)
        frame.pack(fill="both", expand=True, padx=5, pady=5)
        ttk.Label(frame, style="Muted.TLabel",
                  text="📊 Clique em um gráfico para ampliá-lo • ESC ou duplo clique para voltar",
                  ).pack(fill="x", pady=(0, 5))

        self.figure = Figure(figsize=(10, 7), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.figure, master=frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        # Evento do matplotlib (traz .inaxes); o evento Tk não tem esse campo.
        self.canvas.mpl_connect("button_press_event", self._on_graph_click)
        self.plot_placeholder()

    def _build_equations(self):
        frame = ttk.Frame(self.tab_equations, style="Panel.TFrame", padding=22)
        frame.pack(fill="both", expand=True, padx=7, pady=7)
        box = make_text(frame, font=("Consolas", 11))
        sb = ttk.Scrollbar(frame, command=box.yview)
        box.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        box.pack(fill="both", expand=True)
        set_text(box, (
            "EQUAÇÕES UTILIZADAS\n\n"
            "Peso:\nF_g = m · g\n\n"
            "Arrasto aerodinâmico (corpo):\nF_d = 1/2 · ρ · C_d · A · v²\n\n"
            "Arrasto combinado com paraquedas aberto:\n"
            "F_d = 1/2 · ρ · (C_d_corpo·A_corpo + C_d_par·A_par) · v²\n"
            "(corpo e paraquedas contribuem simultaneamente)\n\n"
            "Equação da queda:\nm · dv/dt = m · g - F_d\n"
            "Integração: Runge-Kutta 4ª ordem, dt = 2 ms, com interpolação\n"
            "linear do instante de impacto no último passo.\n\n"
            "Velocidade terminal aproximada:\nv_t = √(2 · m · g / (ρ · C_d · A))\n\n"
            "Energia cinética:\nE_k = 1/2 · m · v²\n\n"
            "Paraquedas (dimensionado só pelo paraquedas, sem o corpo):\n"
            "A_p = 2 · m · g / (ρ · C_d · v_alvo²)\n"
            "→ a velocidade simulada fica um pouco MENOR que a velocidade alvo,\n"
            "  pois o arrasto do corpo também atua.\n\n"
            "Paraquedas circular:  D = 2 · √(A_p / π)\n"
            "Paraquedas retangular:  A = largura · altura\n"
            "Paraquedas cruciforme (modelo simplificado):\n"
            "A ≈ 0,60 · envergadura · altura\n"
            "Retangular e cruciforme usam apenas a PROPORÇÃO largura:altura.\n\n"
            "Abertura do paraquedas: instantânea na altura de abertura (não\n"
            "modela tempo de inflagem nem choque de abertura).\n\n"
            "As cotas e pontos de fixação do desenho são geométricos e devem\n"
            "ser validados fisicamente antes da construção.\n\n"
            "Desaceleração cinemática líquida durante a parada:\n"
            "a_liq = v² / (2 · d)\n\n"
            "Reação estrutural:  N = m · (a_liq + g)\n"
            "Carga em g:  g_load = N / (m · g) = a_liq/g + 1\n\n"
            "Pressão aproximada:  P = N / A\n\n"
            "Energia absorvida por uma camada simplificada:\n"
            "E_abs ≈ σ_crush · A · d_compressão\n\n"
            "NOTA:\nO índice heurístico é um indicador de engenharia baseado nos\n"
            "parâmetros fornecidos. Não é uma probabilidade de sobrevivência."
        ))

    def _build_save(self):
        t = self.tab_save
        t.columnconfigure(0, weight=1)
        t.columnconfigure(1, weight=1)

        fc = self.card(t, "💾 Arquivo do projeto", 0, 0)
        for text, cmd, style in [
            ("💾 Salvar projeto JSON", self.save_project, "Accent.TButton"),
            ("📂 Carregar projeto JSON", self.load_project, "Secondary.TButton"),
            ("📄 Exportar resultados CSV", self.export_csv, "Secondary.TButton"),
        ]:
            ttk.Button(fc, text=text, command=cmd, style=style).pack(fill="x", pady=5)

        oc = self.card(t, "🔎 Otimização", 0, 1)
        self.add_field(oc, "opt_min_padding_mm", "Proteção mínima", "mm")
        self.add_field(oc, "opt_max_padding_mm", "Proteção máxima", "mm")
        self.add_field(oc, "opt_step_mm", "Passo", "mm")
        self.add_field(oc, "opt_max_structure_g", "Massa máxima", "g")

        rc = self.card(t, "🏆 Melhores configurações (menor score = melhor)", 1, 0, 2)
        cols = ("rank", "padding", "mass", "impact", "g", "margin", "geometry")
        heads = ["#", "Proteção", "Massa estrutura", "Impacto", "Carga g",
                 "Margem energ.", "Geometria"]
        self.opt_tree = ttk.Treeview(rc, columns=cols, show="headings", height=12)
        for c, h in zip(cols, heads):
            self.opt_tree.heading(c, text=h)
            self.opt_tree.column(c, width=50 if c == "rank" else 130, anchor="center")
        self.opt_tree.tag_configure("ok", foreground=GREEN)
        self.opt_tree.tag_configure("bad", foreground=RED)
        self.opt_tree.pack(fill="both", expand=True)

        ttk.Button(t, text="🔎 EXECUTAR OTIMIZAÇÃO", style="Accent.TButton",
                   command=self.optimize).grid(
            row=2, column=0, columnspan=2, sticky="ew", padx=7, pady=7)

    # ------------------------------------------------------------------
    # Defaults / modo
    # ------------------------------------------------------------------
    def _set_defaults(self):
        for key, value in DEFAULTS.items():
            if key in self.vars:
                self.vars[key].set(str(value))
        self._refresh_mode()
        self._refresh_parachute_fields()

    def _on_mode_change(self):
        new = self.mode.get()
        if new == "Básico" and self._last_mode == "Avançado":
            changed = [k for k in ADVANCED_KEYS
                       if safe_float(self.vars[k].get()) != safe_float(DEFAULTS[k])]
            if changed and not messagebox.askyesno(
                "Voltar ao modo básico",
                "O modo básico trava os parâmetros avançados nos valores padrão.\n"
                "Os valores avançados alterados serão restaurados. Continuar?",
            ):
                self.mode.set("Avançado")
                return
            for k in ADVANCED_KEYS:
                self.vars[k].set(str(DEFAULTS[k]))
        self._last_mode = new
        self._refresh_mode()

    def _refresh_mode(self):
        advanced = self.mode.get() == "Avançado"
        for key in ADVANCED_KEYS:
            w = self.widgets.get(key)
            if w is not None:
                w.configure(state="normal" if advanced else "disabled")
                self.labels[key].configure(
                    text=self.label_text[key] + ("" if advanced else "  🔒"))
        if advanced:
            self.home_status.configure(
                text="Modo avançado: todos os parâmetros físicos estão disponíveis.")
        else:
            self.home_status.configure(
                text="Modo básico: parâmetros marcados com 🔒 (gravidade, densidade do ar, "
                     "Cd, fator de segurança, Cd do paraquedas, parada estrutural e tensão "
                     "de esmagamento) ficam nos valores padrão.")

    def _on_parachute_type_change(self):
        self._refresh_parachute_fields()
        self.refresh_parachute_drawing()

    def _refresh_parachute_fields(self):
        state = "disabled" if self.vars["parachute_type"].get() == "Redondo" else "normal"
        for key in ("parachute_width_cm", "parachute_height_cm"):
            self.widgets[key].configure(state=state)

    # ------------------------------------------------------------------
    # Coleta e validação
    # ------------------------------------------------------------------
    def num(self, key, gt=None, ge=None, le=None):
        """Lê um campo numérico. Texto inválido gera erro (a 3.3 usava um
        valor padrão silenciosamente)."""
        raw = self.vars[key].get()
        label = self.label_text[key]
        v = safe_float(raw)
        if v is None or not math.isfinite(v):
            raise ValueError(f"'{label}': valor inválido ({raw!r}).")
        if gt is not None and v <= gt:
            raise ValueError(f"'{label}' precisa ser maior que {gt:g}.")
        if ge is not None and v < ge:
            raise ValueError(f"'{label}' não pode ser menor que {ge:g}.")
        if le is not None and v > le:
            raise ValueError(f"'{label}' não pode ser maior que {le:g}.")
        return v

    def collect_data(self):
        n = self.num
        d = {
            "height": n("height_m", gt=0),
            "gravity": n("gravity", gt=0),
            "rho": n("rho", gt=0),
            "wind_speed": n("wind_speed_ms", ge=0),
            "wind_uncertainty": n("wind_uncertainty_pct", ge=0) / 100,
            "egg_mass": n("egg_mass_g", gt=0) / 1000,
            "structure_mass": n("structure_mass_g", ge=0) / 1000,
            "egg_length": n("egg_length_mm", gt=0) / 1000,
            "egg_width": n("egg_width_mm", gt=0) / 1000,
            "egg_height": n("egg_height_mm", gt=0) / 1000,
            "frontal_area": n("frontal_area_cm2", gt=0) / 10000,
            "shell_strength": n("shell_strength_kpa", gt=0) * 1000,
            "allowed_g": n("allowed_g", gt=0),
            "container_length": n("container_length_cm", gt=0) / 100,
            "container_width": n("container_width_cm", gt=0) / 100,
            "container_height": n("container_height_cm", gt=0) / 100,
            "cd": n("cd", gt=0),
            "safety_factor": n("safety_factor", gt=0),
            "para_enabled": self.vars["para_enabled"].get() == "Sim",
            "parachute_type": self.vars["parachute_type"].get(),
            "target_velocity": n("target_velocity_ms", gt=0),
            "deploy_height": n("deploy_height_m", ge=0),
            "para_cd": n("para_cd", gt=0),
            "max_diameter": n("max_diameter_cm", gt=0) / 100,
            "parachute_width": n("parachute_width_cm", gt=0) / 100,
            "parachute_height": n("parachute_height_cm", gt=0) / 100,
            "cord_length": n("cord_length_cm", gt=0) / 100,
            "cord_count": int(n("cord_count", gt=0)),
            "padding": n("padding_mm", ge=0) / 1000,
            "compression": n("compression", ge=0, le=100) / 100,
            "structural_stop": n("structural_stop_mm", ge=0) / 1000,
            "crush_stress": n("crush_stress_kpa", ge=0) * 1000,
            "protection_area": n("protection_area_cm2", gt=0) / 10000,
            "target_margin": n("target_margin", gt=0),
        }

        if d["para_enabled"]:
            if d["deploy_height"] >= d["height"]:
                raise ValueError(
                    "A altura de abertura do paraquedas deve ser menor que a altura da queda.")
            if d["deploy_height"] <= 0:
                raise ValueError(
                    "A altura de abertura deve ser maior que zero quando o paraquedas "
                    "está habilitado (com 0 ele nunca abre antes do impacto).")

        d["total_mass"] = d["egg_mass"] + d["structure_mass"]
        return d

    # ------------------------------------------------------------------
    # Paraquedas
    # ------------------------------------------------------------------
    def _para_geometry(self, d):
        area = parachute_area_for_velocity(
            d["total_mass"], d["gravity"], d["rho"], d["para_cd"], d["target_velocity"])
        ratio = d["parachute_width"] / d["parachute_height"]
        w, h = parachute_dimensions(d["parachute_type"], area, ratio)
        return area, w, h

    def refresh_parachute_drawing(self):
        try:
            d = self.collect_data()
            _, w, h = self._para_geometry(d)
            args = (d["parachute_type"], w, h, d["cord_count"], d["cord_length"])
        except Exception as exc:
            args = ("Redondo", 0.5, 0.5, 4, 0.45)
            self.set_status(f"Desenho do paraquedas com valores padrão: {exc}", YELLOW)
        draw_parachute(self.para_figure, *args)
        self.para_canvas.draw_idle()

    def calculate_parachute(self):
        try:
            d = self.collect_data()
            area, w, h = self._para_geometry(d)
        except Exception as exc:
            messagebox.showerror("Erro no paraquedas", str(exc))
            return

        ptype = d["parachute_type"]
        ratio = d["parachute_width"] / d["parachute_height"]
        feasible = max(w, h) <= d["max_diameter"]
        label = {"Redondo": "Diâmetro", "Retangular": "Largura × altura"}.get(
            ptype, "Envergadura × altura")
        dim_text = (f"{w*100:.2f} cm" if ptype == "Redondo"
                    else f"{w*100:.2f} × {h*100:.2f} cm")

        vt_body = math.sqrt(2 * d["total_mass"] * d["gravity"]
                            / (d["rho"] * d["cd"] * d["frontal_area"]))
        vt_final = math.sqrt(2 * d["total_mass"] * d["gravity"] / (
            d["rho"] * (d["cd"] * d["frontal_area"] + d["para_cd"] * area)))

        set_text(self.para_result, (
            "☂️ PARAQUEDAS CALCULADO\n\n"
            f"Tipo:                       {ptype}\n"
            f"Massa total:                {fmt(d['total_mass']*1000)} g\n"
            f"Velocidade alvo:            {fmt(d['target_velocity'])} m/s\n"
            f"Área necessária:            {fmt(area, 4)} m²\n"
            f"{label}:  {dim_text}\n"
            f"Altura de abertura:         {fmt(d['deploy_height'])} m\n"
            f"Comprimento das cordas:     {fmt(d['cord_length']*100)} cm\n"
            f"Pontos de fixação:          {d['cord_count']}\n"
            f"Velocidade sem paraquedas:  {fmt(vt_body)} m/s\n"
            f"Velocidade terminal c/ par.: {fmt(vt_final)} m/s (corpo + paraquedas)\n\n"
            f"Dimensão máxima informada:  {fmt(d['max_diameter']*100)} cm\n"
            f"Status geométrico: {'✓ CABE' if feasible else '⚠ NÃO CABE'}\n\n"
            "Desenho: vista superior com cotas e pontos de fixação; vista lateral\n"
            "em escala com o comprimento das cordas."
        ))

        for item in self.para_tree.get_children():
            self.para_tree.delete(item)
        for dim_cm in np.arange(20, max(20, d["max_diameter"] * 100) + 0.1, 5):
            a = parachute_area_from_size(ptype, dim_cm / 100, ratio)
            v = math.sqrt(2 * d["total_mass"] * d["gravity"] / (d["rho"] * d["para_cd"] * a))
            e = 0.5 * d["total_mass"] * v * v
            self.para_tree.insert(
                "", "end", tags=("ok" if v <= d["target_velocity"] else "",),
                values=(f"{dim_cm:.0f} cm", f"{a:.3f} m²", f"{v:.2f} m/s", f"{e:.2f} J"))

        self.vars["para_enabled"].set("Sim")
        draw_parachute(self.para_figure, ptype, w, h, d["cord_count"], d["cord_length"])
        self.para_canvas.draw_idle()

        if not feasible:
            messagebox.showwarning(
                "Paraquedas fora da dimensão máxima",
                "A geometria calculada excede a dimensão máxima informada. "
                "Ajuste a dimensão máxima ou a velocidade-alvo.")

        self.run_simulation(navigate=False)
        self.goto("para")   # a 3.3 saltava para a aba Simulação e escondia o resultado

    def _show_fullscreen(self, title, draw_fn):
        win = tk.Toplevel(self.root)
        win.title(title)
        win.configure(bg=BG)
        win.attributes("-fullscreen", True)

        head = tk.Frame(win, bg=BG)
        head.pack(fill="x", padx=18, pady=12)
        tk.Label(head, text=title, bg=BG, fg=TEXT,
                 font=("Segoe UI", 18, "bold")).pack(side="left")
        tk.Label(head, text="ESC ou duplo clique para voltar", bg=BG, fg=MUTED,
                 font=("Segoe UI", 10)).pack(side="right")

        fig = Figure(figsize=(16, 9), dpi=100)
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=18, pady=(0, 18))
        draw_fn(fig)
        canvas.draw()

        close = lambda event=None: win.destroy()
        win.bind("<Escape>", close)
        canvas.get_tk_widget().bind("<Double-Button-1>", close)
        win.focus_force()

    def _open_parachute_fullscreen(self):
        try:
            d = self.collect_data()
            _, w, h = self._para_geometry(d)
            args = (d["parachute_type"], w, h, d["cord_count"], d["cord_length"])
        except Exception as exc:
            messagebox.showwarning("Paraquedas", f"Parâmetros inválidos:\n{exc}")
            return
        self._show_fullscreen("☂️ PROJETO DO PARAQUEDAS",
                              lambda fig: draw_parachute(fig, *args))

    # ------------------------------------------------------------------
    # Simulação / proteção
    # ------------------------------------------------------------------
    def _apply_result(self, result):
        self.result = result
        self.update_simulation_text()
        self.plot_results()

    def run_simulation(self, navigate=True):
        try:
            result = compute_result(self.collect_data())
        except Exception as exc:
            messagebox.showerror("Erro na simulação", str(exc))
            return

        if not result["simulation"]["landed"]:
            messagebox.showwarning(
                "Simulação não convergiu",
                f"O ovo não atingiu o solo em {MAX_SIM_TIME:.0f} s de simulação. "
                "A 'velocidade de impacto' não corresponde a um impacto real — "
                "revise altura, massa, área frontal, gravidade ou paraquedas.")

        self._apply_result(result)
        if navigate:
            self.goto("sim")

    def calculate_protection(self):
        try:
            d = self.collect_data()
            # Sempre recalcula: a 3.3 reaproveitava o último resultado, que
            # podia ser de parâmetros antigos (e usava queda no vácuo sem ele).
            result = compute_result(d)
            capacity, crush_force, comp_dist = energy_absorption_capacity(
                d["padding"], d["compression"], d["crush_stress"], d["protection_area"])
        except Exception as exc:
            messagebox.showerror("Erro no dimensionamento", str(exc))
            return

        self._apply_result(result)
        impact_energy = result["impact"]["impact_energy"]
        margin = capacity / max(impact_energy, 1e-9)

        if margin >= d["target_margin"]:
            verdict = "✓ A capacidade estimada atinge a margem desejada."
        elif margin >= 1:
            verdict = "⚠ Absorve a energia, mas abaixo da margem desejada."
        else:
            verdict = "⚠ A proteção estimada não absorve toda a energia do impacto."

        set_text(self.protection_result, (
            "🧱 DIMENSIONAMENTO\n\n"
            f"Espessura:                 {fmt(d['padding']*1000)} mm\n"
            f"Compressão utilizável:     {fmt(d['compression']*100)} %\n"
            f"Distância de compressão:   {fmt(comp_dist*1000)} mm\n"
            f"Força de esmagamento:      {fmt_si(crush_force)} N\n"
            f"Capacidade de absorção:    {fmt(capacity)} J\n"
            f"Energia a absorver:        {fmt(impact_energy)} J "
            f"(da simulação{' com paraquedas' if d['para_enabled'] else ''})\n"
            f"Margem energética:         {fmt(margin)} x   (desejada: {fmt(d['target_margin'])} x)\n\n"
            f"{verdict}\n\nModelo simplificado: E_abs ≈ σ_crush · A · d."
        ))
        self.goto("protection")

    def update_simulation_text(self):
        if not self.result:
            return
        r = self.result
        d, sim, impact = r["data"], r["simulation"], r["impact"]
        capacity, margin = r["absorption_capacity"], r["energy_margin"]
        target = d["target_margin"]

        checks = {
            "geometria": r["geometry_ok"],
            "energia": margin >= target,
            "força": impact["force_margin"] >= 1,
            "aceleração": impact["g_margin"] >= 1,
        }
        issues = [k for k, ok in checks.items() if not ok]
        idx = impact["heuristic_index"]

        if not sim["landed"]:
            msg, color = "⚠ SIMULAÇÃO NÃO CONVERGIU — resultado não confiável", RED
        elif not issues and idx < 60:
            msg, color = f"✓ APROVADO — risco {impact['risk_label'].lower()}", GREEN
        elif not issues:
            # Dentro de todos os limites, mas perto deles (o índice sozinho
            # contradizia o "aprovado" na 3.3).
            msg, color = (f"✓ DENTRO DOS LIMITES — margem estreita "
                          f"(índice {fmt(idx, 0)}/100)"), YELLOW
        else:
            physical = [i for i in issues if i != "geometria"]
            msg = f"⚠ REPROVADO em: {', '.join(issues)} — risco {impact['risk_label'].lower()}"
            color = RED if physical else YELLOW
        self.sim_status_var.set(msg)
        self.sim_status_label.configure(fg=color)

        self.kpis["v"].configure(text=f"{fmt(sim['impact_velocity'])} m/s")
        self.kpis["g"].configure(
            text=f"{fmt(impact['g_load'], 1)} g",
            fg=GREEN if checks["aceleração"] else RED)
        self.kpis["margin"].configure(
            text=f"{fmt(margin)} x",
            fg=GREEN if checks["energia"] else RED)
        self.kpis["risk"].configure(text=f"{fmt(idx, 0)}/100", fg=color)

        ok = lambda b, yes="✓ OK", no="⚠ FALHA": yes if b else no
        warn = ("" if sim["landed"] else
                "⚠ ATENÇÃO: o ovo não atingiu o solo no tempo máximo. Os valores "
                "de impacto NÃO representam um pouso real.\n\n")
        vt_final = sim["terminal_velocity_final"]

        set_text(self.result_text, (
            warn + "🥚 RESULTADO DA SIMULAÇÃO\n\n"
            f"Massa total:                  {fmt(d['total_mass']*1000)} g\n"
            f"Altura da queda:              {fmt(d['height'])} m\n"
            f"Paraquedas:                   {'SIM' if d['para_enabled'] else 'NÃO'}\n\n"
            "QUEDA\n"
            f"Tempo de queda:               {fmt(sim['fall_time'])} s\n"
            f"Velocidade de impacto:        {fmt(sim['impact_velocity'])} m/s "
            f"({fmt(sim['impact_velocity']*3.6)} km/h)\n"
            f"Velocidade terminal (corpo):  {fmt(sim['terminal_velocity'])} m/s\n"
            + (f"Velocidade terminal (c/ par.): {fmt(vt_final)} m/s\n" if vt_final else "")
            + f"Energia cinética:             {fmt(impact['impact_energy'])} J\n\n"
            "IMPACTO\n"
            f"Distância de parada:          {fmt(impact['stop_distance']*1000)} mm\n"
            f"Desaceleração líquida:        {fmt(impact['net_deceleration'])} m/s²\n"
            f"Carga aproximada (c/ peso):   {fmt(impact['g_load'])} g\n"
            f"Força média na estrutura:     {fmt(impact['force'])} N\n"
            f"Pressão aproximada:           {fmt(impact['pressure']/1000)} kPa\n\n"
            "PROTEÇÃO\n"
            f"Capacidade de absorção:       {fmt(capacity)} J\n"
            f"Margem energética:            {fmt(margin)} x (desejada: {fmt(target)} x)\n"
            f"Geometria:                    {ok(checks['geometria'], '✓ OK', '⚠ NÃO CABE')}\n"
            f"Energia:                      {ok(checks['energia'], '✓ OK', '⚠ INSUFICIENTE')}\n"
            f"Força:                        {ok(checks['força'], '✓ DENTRO', '⚠ EXCEDE')}\n"
            f"Aceleração:                   {ok(checks['aceleração'], '✓ DENTRO', '⚠ EXCEDE')}\n\n"
            f"Índice heurístico:            {fmt(idx)}/100 — {impact['risk_label']}\n\n"
            "⚠ IMPORTANTE\nO índice é heurístico, não uma probabilidade de o ovo "
            "sobreviver. A carga em g já inclui o peso durante a parada "
            "(g_load = a_líquida/g + 1)."
        ))

        summary = (f"Última simulação: {fmt(sim['fall_time'])} s | impacto "
                   f"{fmt(sim['impact_velocity'])} m/s | {msg}")
        self.home_status.configure(text=summary)
        self.set_status(summary, color)

    # ------------------------------------------------------------------
    # Otimização
    # ------------------------------------------------------------------
    def optimize(self):
        try:
            d = self.collect_data()
            min_p = self.num("opt_min_padding_mm", gt=0)
            max_p = self.num("opt_max_padding_mm", gt=0)
            step = self.num("opt_step_mm", gt=0)
            max_mass = self.num("opt_max_structure_g", ge=0)
            if max_p < min_p:
                raise ValueError("Proteção máxima não pode ser menor que a mínima.")
        except Exception as exc:
            messagebox.showerror("Erro na otimização", str(exc))
            return

        self.root.configure(cursor="watch")
        self.root.update_idletasks()
        rows = []
        try:
            for padding_mm in np.arange(min_p, max_p + step / 2, step):
                for structure_g in np.linspace(0, max_mass, 7):
                    d2 = dict(d)
                    d2["padding"] = padding_mm / 1000
                    d2["structure_mass"] = structure_g / 1000
                    d2["total_mass"] = d2["egg_mass"] + d2["structure_mass"]
                    r = compute_result(d2)
                    sim, impact, margin = r["simulation"], r["impact"], r["energy_margin"]

                    score = (impact["heuristic_index"]
                             + max(0, 1 - min(margin, 1)) * 60
                             + (50 if not r["geometry_ok"] else 0)
                             + (1000 if not sim["landed"] else 0))
                    rows.append((score, padding_mm, structure_g, sim["impact_velocity"],
                                 impact["g_load"], margin, r["geometry_ok"], sim["landed"]))
        except Exception as exc:
            messagebox.showerror("Erro na otimização", str(exc))
            return
        finally:
            self.root.configure(cursor="")

        rows.sort(key=lambda x: x[0])
        for item in self.opt_tree.get_children():
            self.opt_tree.delete(item)
        for i, (_, pad, mass, vel, g, margin, geo, landed) in enumerate(rows[:25], 1):
            label = "OK" if geo else "NÃO CABE"
            if not landed:
                label += " (NÃO CONVERGIU)"
            good = geo and landed and margin >= d["target_margin"]
            self.opt_tree.insert(
                "", "end", tags=("ok" if good else "bad" if not geo else "",),
                values=(i, f"{pad:.0f} mm", f"{mass:.1f} g", f"{vel:.2f} m/s",
                        f"{g:.1f} g", f"{margin:.2f} x", label))
        self.goto("save")

    # ------------------------------------------------------------------
    # Gráficos
    # ------------------------------------------------------------------
    def plot_placeholder(self):
        self.figure.clear()
        self.figure.patch.set_facecolor(BG)
        ax = self.figure.add_subplot(111)
        ax.set_facecolor(PANEL)
        ax.text(0.5, 0.5, "Execute uma simulação para gerar os gráficos.",
                ha="center", va="center", color=MUTED, transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_color(BORDER)
        BRAND.watermark_figure(self.figure)
        self.canvas.draw()

    def plot_results(self):
        if not self.result:
            return
        draw_all_graphs(self.figure, self.result["simulation"], self.result["landing"])
        self.canvas.draw()

    def _on_graph_click(self, event):
        if not self.result or event.inaxes is None:
            return
        axes = list(self.figure.axes)
        if event.inaxes not in axes:
            return
        index = axes.index(event.inaxes)
        sim, landing = self.result["simulation"], self.result["landing"]
        self._show_fullscreen(
            "📊 GRÁFICO AMPLIADO",
            lambda fig: draw_single_graph(fig, index, sim, landing))

    # ------------------------------------------------------------------
    # Arquivos
    # ------------------------------------------------------------------
    def project_dict(self):
        return {
            "version": "3.4",
            "mode": self.mode.get(),
            "variables": {k: v.get() for k, v in self.vars.items()},
        }

    def save_project(self):
        path = filedialog.asksaveasfilename(
            title="Salvar projeto", defaultextension=".json",
            filetypes=[("Projeto Egg Drop", "*.json"), ("Todos", "*.*")])
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(self.project_dict(), f, ensure_ascii=False, indent=2)
        except OSError as exc:
            messagebox.showerror("Erro ao salvar", str(exc))
            return
        self.last_project_path = path
        self.set_status(f"Projeto salvo em {path}", GREEN)

    def load_project(self):
        path = filedialog.askopenfilename(
            title="Carregar projeto",
            filetypes=[("Projeto Egg Drop", "*.json"), ("Todos", "*.*")])
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as f:
                project = json.load(f)
            variables = project.get("variables", {})
            if not isinstance(variables, dict):
                raise ValueError("Arquivo de projeto inválido: 'variables' não é um objeto.")
        except (OSError, json.JSONDecodeError, ValueError, AttributeError) as exc:
            messagebox.showerror("Erro ao carregar", str(exc))
            return

        ignored = []
        for key, value in variables.items():
            widget = self.widgets.get(key)
            if widget is None:
                continue
            allowed = widget.cget("values") if isinstance(widget, ttk.Combobox) else None
            if allowed and str(value) not in [str(a) for a in allowed]:
                ignored.append(key)       # valor inválido para um campo de lista
                continue
            self.vars[key].set(str(value))

        if project.get("mode") in ("Básico", "Avançado"):
            self.mode.set(project["mode"])
            self._last_mode = project["mode"]
        self._refresh_mode()
        self._refresh_parachute_fields()
        self.refresh_parachute_drawing()

        # O resultado anterior pertence a outros parâmetros.
        self.result = None
        self.plot_placeholder()
        self.sim_status_var.set("Projeto carregado — execute a simulação.")
        self.sim_status_label.configure(fg=MUTED)
        set_text(self.result_text, "")
        for k in self.kpis.values():
            k.configure(text="—", fg=TEXT)

        self.last_project_path = path
        extra = f" (campos ignorados: {', '.join(ignored)})" if ignored else ""
        self.set_status(f"Projeto carregado de {path}{extra}", GREEN)

    def export_csv(self):
        if not self.result:
            messagebox.showwarning("Nada para exportar", "Execute uma simulação antes de exportar.")
            return
        path = filedialog.asksaveasfilename(
            title="Exportar resultados", defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("Todos", "*.*")])
        if not path:
            return

        r = self.result
        d, sim, impact = r["data"], r["simulation"], r["impact"]
        rows = [
            ("parâmetro", "valor", "unidade"),
            ("massa_total", d["total_mass"], "kg"),
            ("altura", d["height"], "m"),
            ("tempo_queda", sim["fall_time"], "s"),
            ("simulacao_convergiu", sim["landed"], ""),
            ("velocidade_impacto", sim["impact_velocity"], "m/s"),
            ("velocidade_impacto_kmh", sim["impact_velocity"] * 3.6, "km/h"),
            ("energia_impacto", impact["impact_energy"], "J"),
            ("distancia_parada", impact["stop_distance"], "m"),
            ("desaceleracao_liquida", impact["net_deceleration"], "m/s²"),
            ("carga_g", impact["g_load"], "g"),
            ("forca_media", impact["force"], "N"),
            ("pressao", impact["pressure"], "Pa"),
            ("margem_forca", impact["force_margin"], "x"),
            ("margem_g", impact["g_margin"], "x"),
            ("capacidade_absorcao", r["absorption_capacity"], "J"),
            ("margem_energetica", r["energy_margin"], "x"),
            ("indice_heuristico", impact["heuristic_index"], "/100"),
            ("classificacao", impact["risk_label"], ""),
            ("deriva_horizontal_estimada", r["landing"]["position_m"], "m"),
        ]
        # Vírgula decimal: o Excel em pt-BR trata "1.5" como texto/data.
        rows = [tuple(str(c).replace(".", ",") if isinstance(c, float) else c for c in row)
                for row in rows]
        try:
            with open(path, "w", newline="", encoding="utf-8-sig") as f:
                csv.writer(f, delimiter=";").writerows(rows)
        except OSError as exc:
            messagebox.showerror("Erro ao exportar", str(exc))
            return
        self.set_status(f"Resultados exportados para {path}", GREEN)


def main():
    root = tk.Tk()
    EggDropLab(root)
    root.mainloop()


if __name__ == "__main__":
    main()
