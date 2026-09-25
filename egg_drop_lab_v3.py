import csv
import json
import math
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


APP_TITLE = "Egg Drop Lab 3.1"

BG = "#080b12"
PANEL = "#111827"
PANEL_2 = "#182235"
TEXT = "#e5e7eb"
MUTED = "#94a3b8"
PURPLE = "#a78bfa"
GREEN = "#34d399"
YELLOW = "#fbbf24"
RED = "#fb7185"
CYAN = "#67e8f9"

                                                                       
                                                                          
                                                
ADVANCED_KEYS = (
    "gravity",
    "rho",
    "cd",
    "safety_factor",
    "para_cd",
    "structural_stop_mm",
    "crush_stress_kpa",
)


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
    max_time=120.0,
):
    """
    Numerical RK4 simulation.
    y = height above ground, v = downward velocity.
    """
    if height <= 0 or mass <= 0 or area <= 0 or gravity <= 0:
        raise ValueError("Altura, massa, área e gravidade devem ser positivas.")

    t = 0.0
    y = height
    v = 0.0

    ts = [0.0]
    ys = [y]
    vs = [v]
    accs = [gravity]
    energies = [0.0]
    parachute_opened = False
    opening_time = None

    def acceleration(y_, v_):
                                             
        drag_coef_area = cd * area

                                                                           
                                                                     
                                                                           
                                                                           
                                                                         
                                                                        
                                                                         
        if (
            parachute_area > 0
            and parachute_cd > 0
            and y_ <= parachute_deploy_height
        ):
            drag_coef_area += parachute_cd * parachute_area

        drag = 0.5 * rho * drag_coef_area * abs(v_) * v_
        return gravity - drag / mass

    while y > 0 and t < max_time:
        if (
            parachute_area > 0
            and parachute_cd > 0
            and y <= parachute_deploy_height
            and not parachute_opened
        ):
            parachute_opened = True
            opening_time = t

        def dy_dt(y_, v_):
            return -v_

        def dv_dt(y_, v_):
            return acceleration(y_, v_)

        k1y = dy_dt(y, v)
        k1v = dv_dt(y, v)

        k2y = dy_dt(y + 0.5 * dt * k1y, v + 0.5 * dt * k1v)
        k2v = dv_dt(y + 0.5 * dt * k1y, v + 0.5 * dt * k1v)

        k3y = dy_dt(y + 0.5 * dt * k2y, v + 0.5 * dt * k2v)
        k3v = dv_dt(y + 0.5 * dt * k2y, v + 0.5 * dt * k2v)

        k4y = dy_dt(y + dt * k3y, v + dt * k3v)
        k4v = dv_dt(y + dt * k3y, v + dt * k3v)

        y += dt * (k1y + 2 * k2y + 2 * k3y + k4y) / 6
        v += dt * (k1v + 2 * k2v + 2 * k3v + k4v) / 6

        if y < 0:
            y = 0

        t += dt
        a = acceleration(y, v)

        ts.append(t)
        ys.append(y)
        vs.append(v)
        accs.append(a)
        energies.append(0.5 * mass * v * v)

                                                                            
                                                                          
                                                                           
                                                                      
                                                                    
    landed = y <= 0

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
        "terminal_velocity": math.sqrt(
            2 * mass * gravity / (rho * cd * area)
        ),
    }


def parachute_area_for_velocity(mass, gravity, rho, cd, target_velocity):
    if min(mass, gravity, rho, cd, target_velocity) <= 0:
        raise ValueError("Todos os parâmetros do paraquedas devem ser positivos.")
    return 2 * mass * gravity / (rho * cd * target_velocity**2)


def circular_diameter(area):
    return 2 * math.sqrt(area / math.pi)


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

                                                                      
                                                
    net_decel = impact_velocity**2 / (2 * stop_distance)

                                                                       
                                                                           
                                                                         
                                                     
                                                          
                                                                           
                                                                      
                                                                         
                                                                           
                                                                     
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


def energy_absorption_capacity(mass, padding, compression, crush_stress, area):
    """
    Simplified engineering estimate:
    F ~ crush_stress * area
    E_abs = F * compression_distance
    """
    compression_distance = max(0.0, padding * clamp(compression, 0, 1))
    force = max(0.0, crush_stress) * max(0.0, area)
    return force * compression_distance, force, compression_distance


class EggDropLab:
    def __init__(self, root):
        self.root = root
        self.root.title(APP_TITLE)
        self.root.geometry("1420x900")
        self.root.minsize(1180, 760)
        self.root.configure(bg=BG)

        self.vars = {}
        self.widgets = {}
        self.result = None
        self.last_project_path = None
        self.mode = tk.StringVar(value="Básico")

        self._configure_style()
        self._build_header()
        self._build_main()
        self._set_defaults()

                                                              
                
                                                              

    def _configure_style(self):
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass

        style.configure(
            ".",
            background=BG,
            foreground=TEXT,
            fieldbackground=PANEL_2,
            font=("Segoe UI", 10),
        )
        style.configure("TFrame", background=BG)
        style.configure("Panel.TFrame", background=PANEL)
        style.configure("TLabel", background=BG, foreground=TEXT)
        style.configure("Muted.TLabel", background=BG, foreground=MUTED)
        style.configure(
            "Title.TLabel",
            background=BG,
            foreground=TEXT,
            font=("Segoe UI", 20, "bold"),
        )
        style.configure(
            "SubTitle.TLabel",
            background=BG,
            foreground=MUTED,
            font=("Segoe UI", 10),
        )
        style.configure(
            "CardTitle.TLabel",
            background=PANEL,
            foreground=PURPLE,
            font=("Segoe UI", 11, "bold"),
        )
        style.configure(
            "Accent.TButton",
            background=PURPLE,
            foreground="#0b0b12",
            borderwidth=0,
            padding=(14, 9),
            font=("Segoe UI", 10, "bold"),
        )
        style.map(
            "Accent.TButton",
            background=[("active", "#c4b5fd"), ("pressed", "#8b5cf6")],
        )
        style.configure(
            "Secondary.TButton",
            background=PANEL_2,
            foreground=TEXT,
            borderwidth=0,
            padding=(12, 8),
        )
        style.map(
            "Secondary.TButton",
            background=[("active", "#26344d")],
        )
        style.configure(
            "TEntry",
            fieldbackground=PANEL_2,
            foreground=TEXT,
            insertcolor=TEXT,
            borderwidth=0,
            padding=7,
        )
        style.map(
            "TEntry",
            fieldbackground=[("disabled", "#0f1622")],
            foreground=[("disabled", MUTED)],
        )
        style.configure(
            "TCombobox",
            fieldbackground=PANEL_2,
            foreground=TEXT,
            background=PANEL_2,
            borderwidth=0,
            padding=5,
        )
        style.map(
            "TCombobox",
            fieldbackground=[("disabled", "#0f1622")],
            foreground=[("disabled", MUTED)],
        )
        style.configure(
            "TNotebook",
            background=BG,
            borderwidth=0,
        )
        style.configure(
            "TNotebook.Tab",
            background=PANEL,
            foreground=MUTED,
            padding=(14, 8),
        )
        style.map(
            "TNotebook.Tab",
            background=[("selected", PANEL_2)],
            foreground=[("selected", TEXT)],
        )

    def _build_header(self):
        header = ttk.Frame(self.root)
        header.pack(fill="x", padx=20, pady=(16, 8))

        left = ttk.Frame(header)
        left.pack(side="left")

        ttk.Label(left, text="🥚 EGG DROP LAB", style="Title.TLabel").pack(
            anchor="w"
        )
        ttk.Label(
            left,
            text="Simulação, paraquedas e dimensionamento de proteção",
            style="SubTitle.TLabel",
        ).pack(anchor="w")

        right = ttk.Frame(header)
        right.pack(side="right")

        ttk.Label(right, text="Modo:", style="Muted.TLabel").pack(
            side="left", padx=(0, 6)
        )
        mode = ttk.Combobox(
            right,
            textvariable=self.mode,
            values=["Básico", "Avançado"],
            state="readonly",
            width=12,
        )
        mode.pack(side="left")
        mode.bind("<<ComboboxSelected>>", lambda _: self._refresh_mode())

    def _build_main(self):
        container = ttk.Frame(self.root)
        container.pack(fill="both", expand=True, padx=20, pady=(0, 20))

        self.notebook = ttk.Notebook(container)
        self.notebook.pack(fill="both", expand=True)

        self.tab_home = ttk.Frame(self.notebook)
        self.tab_project = ttk.Frame(self.notebook)
        self.tab_para = ttk.Frame(self.notebook)
        self.tab_protection = ttk.Frame(self.notebook)
        self.tab_sim = ttk.Frame(self.notebook)
        self.tab_graphs = ttk.Frame(self.notebook)
        self.tab_equations = ttk.Frame(self.notebook)
        self.tab_save = ttk.Frame(self.notebook)

        self.notebook.add(self.tab_home, text="🏠 Início")
        self.notebook.add(self.tab_project, text="🥚 Projeto")
        self.notebook.add(self.tab_para, text="☂️ Paraquedas")
        self.notebook.add(self.tab_protection, text="🧱 Proteção")
        self.notebook.add(self.tab_sim, text="🧪 Simulação")
        self.notebook.add(self.tab_graphs, text="📊 Gráficos")
        self.notebook.add(self.tab_equations, text="🧮 Equações")
        self.notebook.add(self.tab_save, text="💾 Projeto")

        self._build_home()
        self._build_project()
        self._build_para()
        self._build_protection()
        self._build_sim()
        self._build_graphs()
        self._build_equations()
        self._build_save()

    def card(self, parent, title, row, column=0, columnspan=1):
        frame = ttk.Frame(parent, style="Panel.TFrame", padding=16)
        frame.grid(
            row=row,
            column=column,
            columnspan=columnspan,
            sticky="nsew",
            padx=7,
            pady=7,
        )
        ttk.Label(frame, text=title, style="CardTitle.TLabel").pack(
            anchor="w", pady=(0, 10)
        )
        return frame

    def add_field(self, parent, key, label, default="", unit="", values=None):
        line = ttk.Frame(parent, style="Panel.TFrame")
        line.pack(fill="x", pady=4)

        ttk.Label(line, text=label).pack(side="left", fill="x", expand=True)

        var = tk.StringVar(value=str(default))
        self.vars[key] = var

        if values:
            widget = ttk.Combobox(
                line,
                textvariable=var,
                values=values,
                state="readonly",
                width=13,
            )
        else:
            widget = ttk.Entry(line, textvariable=var, width=14)

        widget.pack(side="right")
        if unit:
            ttk.Label(
                line,
                text=unit,
                style="Muted.TLabel",
                width=8,
                anchor="w",
            ).pack(side="right", padx=(7, 0))

                                                                      
                                                                      
        self.widgets[key] = widget

        return widget

    def add_section_button(self, parent, text, command, row, column=0):
        button = ttk.Button(
            parent,
            text=text,
            command=command,
            style="Accent.TButton",
        )
        button.grid(
            row=row,
            column=column,
            sticky="ew",
            padx=7,
            pady=7,
        )
        return button

                                                              
          
                                                              

    def _build_home(self):
        self.tab_home.columnconfigure(0, weight=1)
        self.tab_home.columnconfigure(1, weight=1)

        hero = ttk.Frame(self.tab_home, style="Panel.TFrame", padding=28)
        hero.grid(row=0, column=0, columnspan=2, sticky="ew", padx=7, pady=7)

        ttk.Label(
            hero,
            text="🥚 Egg Drop Lab",
            background=PANEL,
            foreground=TEXT,
            font=("Segoe UI", 28, "bold"),
        ).pack(anchor="w")
        ttk.Label(
            hero,
            text=(
                "Projete, simule e compare soluções para proteger um ovo "
                "durante uma queda."
            ),
            background=PANEL,
            foreground=MUTED,
            font=("Segoe UI", 11),
        ).pack(anchor="w", pady=(4, 15))

        buttons = ttk.Frame(hero, style="Panel.TFrame")
        buttons.pack(fill="x")

        actions = [
            ("🚀 Criar / simular projeto", self.go_simulation),
            ("☂️ Calcular paraquedas", self.go_parachute),
            ("🧱 Dimensionar proteção", self.go_protection),
            ("📊 Analisar projeto", self.go_simulation),
        ]

        for i, (text, command) in enumerate(actions):
            ttk.Button(
                buttons,
                text=text,
                command=command,
                style="Accent.TButton",
            ).grid(row=0, column=i, sticky="ew", padx=4)
            buttons.columnconfigure(i, weight=1)

        info = self.card(
            self.tab_home,
            "Como funciona",
            1,
            0,
            2,
        )
        ttk.Label(
            info,
            text=(
                "1. Informe os dados do ovo e da estrutura.\n"
                "2. Se quiser, ative o paraquedas e peça o dimensionamento automático.\n"
                "3. Simule a queda.\n"
                "4. Veja velocidade, aceleração, energia e margem de segurança.\n\n"
                "⚠️ Os resultados são estimativas de engenharia. "
                "Eles não garantem que um ovo real sobreviverá ao teste."
            ),
            justify="left",
            background=PANEL,
            foreground=TEXT,
        ).pack(anchor="w")

        self.home_status = ttk.Label(
            self.tab_home,
            text="Nenhuma simulação executada.",
            style="Muted.TLabel",
        )
        self.home_status.grid(
            row=2,
            column=0,
            columnspan=2,
            sticky="w",
            padx=14,
            pady=10,
        )

    def _build_project(self):
        self.tab_project.columnconfigure(0, weight=1)
        self.tab_project.columnconfigure(1, weight=1)

        env = self.card(self.tab_project, "🌍 Queda e ambiente", 0, 0)
        self.add_field(env, "height_m", "Altura da queda", 3.0, "m")
        self.add_field(env, "gravity", "Gravidade", 9.80665, "m/s²")
        self.add_field(env, "rho", "Densidade do ar", 1.225, "kg/m³")

        egg = self.card(self.tab_project, "🥚 Ovo", 0, 1)
        self.add_field(egg, "egg_mass_g", "Massa", 60, "g")
        self.add_field(egg, "egg_length_mm", "Comprimento", 60, "mm")
        self.add_field(egg, "egg_width_mm", "Largura", 45, "mm")
        self.add_field(egg, "egg_height_mm", "Altura", 45, "mm")
        self.add_field(egg, "frontal_area_cm2", "Área frontal", 20, "cm²")
        self.add_field(egg, "shell_strength_kpa", "Resistência estimada", 300, "kPa")
        self.add_field(egg, "allowed_g", "Aceleração limite", 100, "g")

        structure = self.card(
            self.tab_project,
            "📦 Estrutura / recipiente",
            1,
            0,
        )
        self.add_field(structure, "structure_mass_g", "Massa da estrutura", 100, "g")
        self.add_field(structure, "container_length_cm", "Comprimento interno", 15, "cm")
        self.add_field(structure, "container_width_cm", "Largura interna", 15, "cm")
        self.add_field(structure, "container_height_cm", "Altura interna", 15, "cm")

        aero = self.card(self.tab_project, "🌬️ Aerodinâmica", 1, 1)
        self.add_field(aero, "cd", "Coeficiente de arrasto", 0.47, "")
        self.add_field(aero, "safety_factor", "Fator de segurança", 2.0, "")

        buttons = ttk.Frame(self.tab_project)
        buttons.grid(row=2, column=0, columnspan=2, sticky="ew", pady=8)
        buttons.columnconfigure(0, weight=1)
        buttons.columnconfigure(1, weight=1)

        ttk.Button(
            buttons,
            text="🧮 Calcular / simular",
            command=self.run_simulation,
            style="Accent.TButton",
        ).grid(row=0, column=0, sticky="ew", padx=7)
        ttk.Button(
            buttons,
            text="☂️ Dimensionar paraquedas",
            command=self.calculate_parachute,
            style="Secondary.TButton",
        ).grid(row=0, column=1, sticky="ew", padx=7)

    def _build_para(self):
        self.tab_para.columnconfigure(0, weight=1)
        self.tab_para.columnconfigure(1, weight=1)

        p = self.card(self.tab_para, "☂️ Calculadora de paraquedas", 0, 0)
        self.add_field(p, "para_enabled", "Usar paraquedas", "Não", "", ["Não", "Sim"])
        self.add_field(p, "target_velocity_ms", "Velocidade alvo", 4.0, "m/s")
        self.add_field(p, "deploy_height_m", "Altura de abertura", 1.5, "m")
        self.add_field(p, "para_cd", "Cd do paraquedas", 1.5, "")
        self.add_field(p, "max_diameter_cm", "Diâmetro máximo", 60, "cm")

        p2 = self.card(self.tab_para, "📐 Resultado", 0, 1)
        self.para_result = tk.Text(
            p2,
            height=13,
            width=55,
            bg=PANEL_2,
            fg=TEXT,
            insertbackground=TEXT,
            relief="flat",
            font=("Consolas", 10),
            padx=12,
            pady=10,
        )
        self.para_result.pack(fill="both", expand=True)

        table_frame = self.card(
            self.tab_para,
            "🔎 Alternativas de diâmetro",
            1,
            0,
            2,
        )

        columns = ("diameter", "area", "velocity", "energy")
        self.para_tree = ttk.Treeview(
            table_frame,
            columns=columns,
            show="headings",
            height=9,
        )
        headers = {
            "diameter": "Diâmetro",
            "area": "Área",
            "velocity": "Velocidade",
            "energy": "Energia",
        }
        for c in columns:
            self.para_tree.heading(c, text=headers[c])
            self.para_tree.column(c, width=170, anchor="center")
        self.para_tree.pack(fill="both", expand=True)

        ttk.Button(
            self.tab_para,
            text="☂️ CALCULAR PARAQUEDAS",
            command=self.calculate_parachute,
            style="Accent.TButton",
        ).grid(row=2, column=0, columnspan=2, sticky="ew", padx=7, pady=8)

    def _build_protection(self):
        self.tab_protection.columnconfigure(0, weight=1)
        self.tab_protection.columnconfigure(1, weight=1)

        prot = self.card(self.tab_protection, "🧱 Amortecimento", 0, 0)
        self.add_field(prot, "padding_mm", "Espessura da proteção", 40, "mm")
        self.add_field(prot, "compression", "Compressão utilizável", 70, "%")
        self.add_field(prot, "structural_stop_mm", "Parada estrutural", 5, "mm")
        self.add_field(prot, "crush_stress_kpa", "Tensão de esmagamento", 80, "kPa")
        self.add_field(prot, "protection_area_cm2", "Área efetiva", 30, "cm²")

        limits = self.card(self.tab_protection, "🎯 Limites", 0, 1)
        self.add_field(limits, "target_margin", "Margem desejada", 1.20, "x")

        self.protection_result = tk.Text(
            limits,
            height=11,
            bg=PANEL_2,
            fg=TEXT,
            insertbackground=TEXT,
            relief="flat",
            font=("Consolas", 10),
            padx=12,
            pady=10,
        )
        self.protection_result.pack(fill="both", expand=True)

        ttk.Button(
            self.tab_protection,
            text="🧱 DIMENSIONAR PROTEÇÃO",
            command=self.calculate_protection,
            style="Accent.TButton",
        ).grid(row=1, column=0, columnspan=2, sticky="ew", padx=7, pady=8)

    def _build_sim(self):
        self.tab_sim.columnconfigure(0, weight=1)
        self.tab_sim.columnconfigure(1, weight=1)

                                                                       
                                                                         
        self.sim_status_var = tk.StringVar(value="Nenhuma simulação executada.")
        self.sim_status_label = tk.Label(
            self.tab_sim,
            textvariable=self.sim_status_var,
            bg=PANEL,
            fg=MUTED,
            font=("Segoe UI", 12, "bold"),
            anchor="w",
            padx=16,
            pady=10,
        )
        self.sim_status_label.grid(
            row=0, column=0, columnspan=2, sticky="ew", padx=7, pady=(7, 0)
        )

        self.sim_summary = self.card(
            self.tab_sim,
            "📋 Resultado da simulação",
            1,
            0,
            2,
        )

        self.result_text = tk.Text(
            self.sim_summary,
            height=19,
            bg=PANEL_2,
            fg=TEXT,
            insertbackground=TEXT,
            relief="flat",
            font=("Consolas", 10),
            padx=14,
            pady=12,
        )
        self.result_text.pack(fill="both", expand=True)

        actions = ttk.Frame(self.tab_sim)
        actions.grid(row=2, column=0, columnspan=2, sticky="ew", pady=7)
        actions.columnconfigure(0, weight=1)
        actions.columnconfigure(1, weight=1)
        actions.columnconfigure(2, weight=1)

        ttk.Button(
            actions,
            text="🧪 SIMULAR",
            command=self.run_simulation,
            style="Accent.TButton",
        ).grid(row=0, column=0, sticky="ew", padx=7)
        ttk.Button(
            actions,
            text="🔎 OTIMIZAR",
            command=self.optimize,
            style="Secondary.TButton",
        ).grid(row=0, column=1, sticky="ew", padx=7)
        ttk.Button(
            actions,
            text="📊 VER GRÁFICOS",
            command=self.go_graphs,
            style="Secondary.TButton",
        ).grid(row=0, column=2, sticky="ew", padx=7)

    def _build_graphs(self):
        self.graph_frame = ttk.Frame(self.tab_graphs)
        self.graph_frame.pack(fill="both", expand=True, padx=5, pady=5)

        self.figure = Figure(figsize=(10, 7), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.graph_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.plot_placeholder()

    def _build_equations(self):
        frame = ttk.Frame(self.tab_equations, style="Panel.TFrame", padding=22)
        frame.pack(fill="both", expand=True, padx=7, pady=7)

        text = (
            "EQUAÇÕES UTILIZADAS\n\n"
            "Peso:\n"
            "F_g = m · g\n\n"
            "Arrasto aerodinâmico (corpo):\n"
            "F_d = 1/2 · ρ · C_d · A · v²\n\n"
            "Arrasto combinado com paraquedas aberto:\n"
            "F_d = 1/2 · ρ · (C_d_corpo·A_corpo + C_d_par·A_par) · v²\n"
            "(corpo e paraquedas contribuem simultaneamente para o arrasto)\n\n"
            "Equação da queda:\n"
            "m · dv/dt = m · g - F_d\n\n"
            "Velocidade terminal aproximada:\n"
            "v_t = √(2 · m · g / (ρ · C_d · A))\n\n"
            "Energia cinética:\n"
            "E_k = 1/2 · m · v²\n\n"
            "Paraquedas:\n"
            "A_p = 2 · m · g / (ρ · C_d · v_alvo²)\n\n"
            "Paraquedas circular:\n"
            "D = 2 · √(A_p / π)\n\n"
            "Desaceleração cinemática líquida durante a parada:\n"
            "a_liq = v² / (2 · d)\n\n"
            "Reação estrutural (força que a proteção/casca sofre):\n"
            "N = m · (a_liq + g)\n"
            "(o peso continua atuando durante a parada; N ≠ m·a_liq)\n\n"
            "Carga em g (convenção de acelerômetro):\n"
            "g_load = N / (m · g) = a_liq/g + 1\n\n"
            "Pressão aproximada:\n"
            "P = N / A\n\n"
            "Energia absorvida por uma camada simplificada:\n"
            "E_abs ≈ σ_crush · A · d_compressão\n\n"
            "NOTA:\n"
            "O índice de risco/margem é um indicador de engenharia baseado "
            "nos parâmetros fornecidos. Não é uma probabilidade estatística "
            "de sobrevivência do ovo."
        )

        box = tk.Text(
            frame,
            bg=PANEL_2,
            fg=TEXT,
            relief="flat",
            font=("Consolas", 11),
            padx=18,
            pady=18,
            wrap="word",
        )
        box.insert("1.0", text)
        box.configure(state="disabled")
        box.pack(fill="both", expand=True)

    def _build_save(self):
        self.tab_save.columnconfigure(0, weight=1)
        self.tab_save.columnconfigure(1, weight=1)

        file_card = self.card(self.tab_save, "💾 Arquivo do projeto", 0, 0)
        ttk.Button(
            file_card,
            text="💾 Salvar projeto JSON",
            command=self.save_project,
            style="Accent.TButton",
        ).pack(fill="x", pady=5)
        ttk.Button(
            file_card,
            text="📂 Carregar projeto JSON",
            command=self.load_project,
            style="Secondary.TButton",
        ).pack(fill="x", pady=5)
        ttk.Button(
            file_card,
            text="📄 Exportar resultados CSV",
            command=self.export_csv,
            style="Secondary.TButton",
        ).pack(fill="x", pady=5)

        opt_card = self.card(self.tab_save, "🔎 Otimização", 0, 1)
        self.add_field(opt_card, "opt_min_padding_mm", "Proteção mínima", 5, "mm")
        self.add_field(opt_card, "opt_max_padding_mm", "Proteção máxima", 80, "mm")
        self.add_field(opt_card, "opt_step_mm", "Passo", 5, "mm")
        self.add_field(opt_card, "opt_max_structure_g", "Massa máxima", 300, "g")

        self.optimization_text = tk.Text(
            self.tab_save,
            height=15,
            bg=PANEL_2,
            fg=TEXT,
            relief="flat",
            font=("Consolas", 10),
        )
        self.optimization_text.grid(
            row=1,
            column=0,
            columnspan=2,
            sticky="nsew",
            padx=7,
            pady=7,
        )
        self.tab_save.rowconfigure(1, weight=1)

        ttk.Button(
            self.tab_save,
            text="🔎 EXECUTAR OTIMIZAÇÃO",
            command=self.optimize,
            style="Accent.TButton",
        ).grid(row=2, column=0, columnspan=2, sticky="ew", padx=7, pady=7)

                                                              
                     
                                                              

    def _set_defaults(self):
        defaults = {
            "height_m": 3.0,
            "gravity": 9.80665,
            "rho": 1.225,
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
            "target_velocity_ms": 4.0,
            "deploy_height_m": 1.5,
            "para_cd": 1.5,
            "max_diameter_cm": 60,
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
        for key, value in defaults.items():
            if key in self.vars:
                self.vars[key].set(str(value))

        self._refresh_mode()

    def _refresh_mode(self):
        advanced = self.mode.get() == "Avançado"

                                                                            
                                                                   
                                                                           
                                                                         
                                                                      
                                           
        for key in ADVANCED_KEYS:
            widget = self.widgets.get(key)
            if widget is not None:
                widget.configure(state=("normal" if advanced else "disabled"))

        if advanced:
            self.home_status.configure(
                text="Modo avançado: todos os parâmetros físicos estão disponíveis."
            )
        else:
            self.home_status.configure(
                text=(
                    "Modo básico: gravidade, densidade do ar, Cd, fator de "
                    "segurança, Cd do paraquedas, parada estrutural e tensão "
                    "de esmagamento ficam travados nos valores padrão."
                )
            )

                                                              
                     
                                                              

    def collect_data(self):
        def get(key, default=None):
            return safe_float(self.vars[key].get(), default)

        height = get("height_m")
        gravity = get("gravity")
        rho = get("rho")
        egg_mass = get("egg_mass_g")
        structure_mass = get("structure_mass_g")

        required = {
            "Altura": height,
            "Gravidade": gravity,
            "Densidade do ar": rho,
            "Massa do ovo": egg_mass,
            "Massa da estrutura": structure_mass,
        }

        for name, value in required.items():
            if value is None or value <= 0:
                raise ValueError(f"{name} precisa ser maior que zero.")

        data = {
            "height": height,
            "gravity": gravity,
            "rho": rho,
            "egg_mass": egg_mass / 1000,
            "structure_mass": structure_mass / 1000,
            "egg_length": get("egg_length_mm", 60) / 1000,
            "egg_width": get("egg_width_mm", 45) / 1000,
            "egg_height": get("egg_height_mm", 45) / 1000,
            "frontal_area": get("frontal_area_cm2", 20) / 10000,
            "shell_strength": get("shell_strength_kpa", 300) * 1000,
            "allowed_g": get("allowed_g", 100),
            "container_length": get("container_length_cm", 15) / 100,
            "container_width": get("container_width_cm", 15) / 100,
            "container_height": get("container_height_cm", 15) / 100,
            "cd": get("cd", 0.47),
            "safety_factor": get("safety_factor", 2.0),
            "para_enabled": self.vars["para_enabled"].get() == "Sim",
            "target_velocity": get("target_velocity_ms", 4),
            "deploy_height": get("deploy_height_m", 1.5),
            "para_cd": get("para_cd", 1.5),
            "max_diameter": get("max_diameter_cm", 60) / 100,
            "padding": get("padding_mm", 40) / 1000,
            "compression": get("compression", 70) / 100,
            "structural_stop": get("structural_stop_mm", 5) / 1000,
            "crush_stress": get("crush_stress_kpa", 80) * 1000,
            "protection_area": get("protection_area_cm2", 30) / 10000,
        }

        if data["frontal_area"] <= 0:
            raise ValueError("A área frontal deve ser maior que zero.")
        if data["cd"] <= 0:
            raise ValueError("O coeficiente de arrasto deve ser maior que zero.")
        if not 0 <= data["compression"] <= 1:
            raise ValueError("A compressão deve estar entre 0% e 100%.")

                                                                         
                                                                          
                                                                      
                                                                     
                                                                 
                                                                
                                    
        if data["para_enabled"] and data["deploy_height"] >= data["height"]:
            raise ValueError(
                "A altura de abertura do paraquedas deve ser menor que a altura da queda."
            )
        if data["para_enabled"] and data["deploy_height"] <= 0:
            raise ValueError(
                "A altura de abertura do paraquedas deve ser maior que zero "
                "quando o paraquedas está habilitado (com 0 ele nunca abre "
                "antes do impacto)."
            )

        data["total_mass"] = data["egg_mass"] + data["structure_mass"]

        return data

                                                              
               
                                                              

    def calculate_parachute(self):
        try:
            d = self.collect_data()

            area = parachute_area_for_velocity(
                d["total_mass"],
                d["gravity"],
                d["rho"],
                d["para_cd"],
                d["target_velocity"],
            )
            diameter = circular_diameter(area)

            vt_without = math.sqrt(
                2 * d["total_mass"] * d["gravity"]
                / (d["rho"] * d["cd"] * d["frontal_area"])
            )

            max_diameter = d["max_diameter"]
            feasible = diameter <= max_diameter

                                                  
            candidates = []
            for diameter_cm in np.arange(20, max(20, max_diameter * 100) + 0.1, 5):
                diameter_m = diameter_cm / 100
                candidate_area = math.pi * (diameter_m / 2) ** 2
                velocity = math.sqrt(
                    2 * d["total_mass"] * d["gravity"]
                    / (d["rho"] * d["para_cd"] * candidate_area)
                )
                energy = 0.5 * d["total_mass"] * velocity**2
                candidates.append(
                    (diameter_cm, candidate_area, velocity, energy)
                )

            self.para_result.delete("1.0", "end")
            self.para_result.insert(
                "end",
                "☂️ PARAQUEDAS CALCULADO\n\n"
                f"Massa total:              {fmt(d['total_mass']*1000)} g\n"
                f"Velocidade alvo:          {fmt(d['target_velocity'])} m/s\n"
                f"Área necessária:          {fmt(area, 4)} m²\n"
                f"Diâmetro circular:        {fmt(diameter*100)} cm\n"
                f"Altura de abertura:       {fmt(d['deploy_height'])} m\n"
                f"Velocidade sem paraquedas: {fmt(vt_without)} m/s\n\n"
                f"Diâmetro máximo informado: {fmt(max_diameter*100)} cm\n"
                f"Status: {'✓ CABE' if feasible else '⚠ NÃO CABE'}\n\n"
                "Sugestão:\n"
                + (
                    "O diâmetro calculado está dentro do limite informado."
                    if feasible
                    else "Aumente o limite de diâmetro ou aceite uma velocidade-alvo maior."
                ),
            )

            for item in self.para_tree.get_children():
                self.para_tree.delete(item)

            for diameter_cm, candidate_area, velocity, energy in candidates:
                tag = "ok" if velocity <= d["target_velocity"] else ""
                self.para_tree.insert(
                    "",
                    "end",
                    values=(
                        f"{diameter_cm:.0f} cm",
                        f"{candidate_area:.3f} m²",
                        f"{velocity:.2f} m/s",
                        f"{energy:.2f} J",
                    ),
                    tags=(tag,),
                )

            self.para_tree.tag_configure("ok", foreground=GREEN)

                                                                            
                                                                    
                                                                        
                                                                          
                                                                         
                                                                           
                                                                      
            self.vars["para_enabled"].set("Sim")
            if not feasible:
                messagebox.showwarning(
                    "Paraquedas não cabe no limite informado",
                    "O diâmetro necessário para atingir a velocidade alvo "
                    f"é {fmt(diameter*100)} cm, mas o diâmetro máximo "
                    f"informado é {fmt(max_diameter*100)} cm.\n\n"
                    "Aumente o limite de diâmetro ou aceite uma velocidade "
                    "de impacto maior antes de prosseguir.",
                )
            self.run_simulation()

        except Exception as exc:
            messagebox.showerror("Erro no paraquedas", str(exc))

                                                              
                
                                                              

    def calculate_protection(self):
        try:
            d = self.collect_data()
            capacity, crush_force, compression_distance = energy_absorption_capacity(
                d["total_mass"],
                d["padding"],
                d["compression"],
                d["crush_stress"],
                d["protection_area"],
            )

                                                                   
            if self.result:
                impact_energy = self.result["impact"]["impact_energy"]
            else:
                impact_energy = 0.5 * d["total_mass"] * (
                    math.sqrt(2 * d["gravity"] * d["height"])
                ) ** 2

            margin = capacity / max(impact_energy, 1e-9)

            self.protection_result.delete("1.0", "end")
            self.protection_result.insert(
                "end",
                "🧱 DIMENSIONAMENTO\n\n"
                f"Espessura:                 {fmt(d['padding']*1000)} mm\n"
                f"Compressão utilizável:     {fmt(d['compression']*100)} %\n"
                f"Distância de compressão:   {fmt(compression_distance*1000)} mm\n"
                f"Força de esmagamento:      {fmt_si(crush_force)} N\n"
                f"Capacidade de absorção:    {fmt(capacity)} J\n"
                f"Energia a absorver:        {fmt(impact_energy)} J\n"
                f"Margem energética:         {fmt(margin, 2)} x\n\n"
                + (
                    "✓ A capacidade estimada supera a energia do impacto."
                    if margin >= 1
                    else "⚠ A proteção estimada não absorve toda a energia do impacto."
                )
                + "\n\n"
                "Modelo simplificado: E_abs ≈ σ_crush · A · d."
            )
            self.notebook.select(self.tab_protection)

        except Exception as exc:
            messagebox.showerror("Erro no dimensionamento", str(exc))

                                                              
                
                                                              

    def run_simulation(self):
        try:
            d = self.collect_data()

            parachute_area = 0.0
            if d["para_enabled"]:
                parachute_area = parachute_area_for_velocity(
                    d["total_mass"],
                    d["gravity"],
                    d["rho"],
                    d["para_cd"],
                    d["target_velocity"],
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

            if not sim["landed"]:
                messagebox.showwarning(
                    "Simulação não convergiu",
                    "O ovo não atingiu o solo dentro do tempo máximo de "
                    "simulação (120 s). A 'velocidade de impacto' exibida "
                    "não corresponde a um impacto real — revise altura, "
                    "massa, área frontal, gravidade ou parâmetros do "
                    "paraquedas.",
                )

            stop_distance = (
                d["padding"] * d["compression"]
                + d["structural_stop"]
            )

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
                d["total_mass"],
                d["padding"],
                d["compression"],
                d["crush_stress"],
                d["protection_area"],
            )

            geometry_ok = (
                d["egg_length"] + 2 * d["padding"]
                <= d["container_length"]
                and d["egg_width"] + 2 * d["padding"]
                <= d["container_width"]
                and d["egg_height"] + 2 * d["padding"]
                <= d["container_height"]
            )

            energy_margin = capacity / max(impact["impact_energy"], 1e-9)

            result = {
                "data": d,
                "simulation": sim,
                "impact": impact,
                "absorption_capacity": capacity,
                "energy_margin": energy_margin,
                "geometry_ok": geometry_ok,
            }

            self.result = result
            self.update_simulation_text()
            self.plot_results()
            self.notebook.select(self.tab_sim)

            status = (
                "✓ margem energética adequada"
                if energy_margin >= 1
                else "⚠ capacidade de absorção insuficiente"
            )

            self.home_status.configure(
                text=(
                    f"Última simulação: {fmt(sim['fall_time'])} s | "
                    f"impacto {fmt(sim['impact_velocity'])} m/s | {status}"
                )
            )

        except Exception as exc:
            messagebox.showerror("Erro na simulação", str(exc))

    def update_simulation_text(self):
        if not self.result:
            return

        d = self.result["data"]
        sim = self.result["simulation"]
        impact = self.result["impact"]
        capacity = self.result["absorption_capacity"]
        margin = self.result["energy_margin"]

        para = d["para_enabled"]

        status_geometry = "✓ OK" if self.result["geometry_ok"] else "⚠ NÃO CABE"
        status_energy = "✓ OK" if margin >= 1 else "⚠ INSUFICIENTE"
        status_force = "✓ DENTRO" if impact["force_margin"] >= 1 else "⚠ EXCEDE"
        status_g = "✓ DENTRO" if impact["g_margin"] >= 1 else "⚠ EXCEDE"

        all_ok = (
            self.result["geometry_ok"]
            and margin >= 1
            and impact["force_margin"] >= 1
            and impact["g_margin"] >= 1
            and sim["landed"]
        )

                                               
        if not sim["landed"]:
            self.sim_status_var.set("⚠ SIMULAÇÃO NÃO CONVERGIU — resultado não confiável")
            self.sim_status_label.configure(fg=RED)
        elif all_ok:
            self.sim_status_var.set(f"✓ APROVADO — risco {impact['risk_label'].lower()}")
            self.sim_status_label.configure(fg=GREEN)
        elif impact["heuristic_index"] >= 60:
            self.sim_status_var.set(f"⚠ RISCO {impact['risk_label']} — revise a proteção")
            self.sim_status_label.configure(fg=RED)
        else:
            self.sim_status_var.set(f"⚠ ATENÇÃO — risco {impact['risk_label'].lower()}")
            self.sim_status_label.configure(fg=YELLOW)

        warning_block = ""
        if not sim["landed"]:
            warning_block = (
                "⚠ ATENÇÃO: o ovo não atingiu o solo dentro do tempo máximo "
                "de simulação. Os valores de impacto abaixo NÃO representam "
                "um pouso real.\n\n"
            )

        text = (
            warning_block
            + "🥚 RESULTADO DA SIMULAÇÃO\n\n"
            f"Massa total:                 {fmt(d['total_mass']*1000)} g\n"
            f"Altura da queda:             {fmt(d['height'])} m\n"
            f"Paraquedas:                  {'SIM' if para else 'NÃO'}\n\n"
            "QUEDA\n"
            f"Tempo de queda:              {fmt(sim['fall_time'])} s\n"
            f"Velocidade de impacto:       {fmt(sim['impact_velocity'])} m/s\n"
            f"Velocidade de impacto:       {fmt(sim['impact_velocity']*3.6)} km/h\n"
            f"Velocidade terminal estimada: {fmt(sim['terminal_velocity'])} m/s\n"
            f"Energia cinética:             {fmt(impact['impact_energy'])} J\n\n"
            "IMPACTO\n"
            f"Distância de parada:          {fmt(impact['stop_distance']*1000)} mm\n"
            f"Desaceleração líquida:        {fmt(impact['net_deceleration'])} m/s²\n"
            f"Carga aproximada (c/ peso):   {fmt(impact['g_load'])} g\n"
            f"Força média na estrutura:     {fmt(impact['force'])} N\n"
            f"Pressão aproximada:           {fmt(impact['pressure']/1000)} kPa\n\n"
            "PROTEÇÃO\n"
            f"Capacidade de absorção:       {fmt(capacity)} J\n"
            f"Margem energética:            {fmt(margin, 2)} x\n"
            f"Geometria:                    {status_geometry}\n"
            f"Energia:                      {status_energy}\n"
            f"Força:                        {status_force}\n"
            f"Aceleração:                   {status_g}\n\n"
            f"Índice heurístico:            {fmt(impact['heuristic_index'])}/100\n"
            f"Classificação:                {impact['risk_label']}\n\n"
            "⚠ IMPORTANTE\n"
            "O índice acima é um indicador heurístico, não uma "
            "probabilidade estatística de o ovo sobreviver. A carga em g "
            "já inclui o efeito do peso durante a parada (convenção de "
            "acelerômetro: g_load = a_líquida/g + 1)."
        )

        self.result_text.delete("1.0", "end")
        self.result_text.insert("end", text)

                                                              
                  
                                                              

    def optimize(self):
        try:
            d = self.collect_data()

            min_p = safe_float(self.vars["opt_min_padding_mm"].get(), 5)
            max_p = safe_float(self.vars["opt_max_padding_mm"].get(), 80)
            step = safe_float(self.vars["opt_step_mm"].get(), 5)
            max_mass = safe_float(
                self.vars["opt_max_structure_g"].get(), 300
            )

            if min_p <= 0 or max_p < min_p or step <= 0:
                raise ValueError("Valores de otimização inválidos.")

            rows = []

            for padding_mm in np.arange(min_p, max_p + step / 2, step):
                for structure_g in np.linspace(
                    0, max_mass, 7
                ):
                    d2 = dict(d)
                    d2["padding"] = padding_mm / 1000
                    d2["structure_mass"] = structure_g / 1000
                    d2["total_mass"] = d2["egg_mass"] + d2["structure_mass"]

                    para_area = 0.0
                    if d2["para_enabled"]:
                        para_area = parachute_area_for_velocity(
                            d2["total_mass"],
                            d2["gravity"],
                            d2["rho"],
                            d2["para_cd"],
                            d2["target_velocity"],
                        )

                    sim = simulate_fall(
                        d2["height"],
                        d2["total_mass"],
                        d2["frontal_area"],
                        d2["cd"],
                        d2["rho"],
                        d2["gravity"],
                        para_area,
                        d2["para_cd"],
                        d2["deploy_height"],
                    )

                    impact = impact_analysis(
                        d2["total_mass"],
                        sim["impact_velocity"],
                        d2["padding"] * d2["compression"]
                        + d2["structural_stop"],
                        d2["frontal_area"],
                        d2["shell_strength"],
                        d2["safety_factor"],
                        d2["allowed_g"],
                        d2["gravity"],
                    )

                    capacity, _, _ = energy_absorption_capacity(
                        d2["total_mass"],
                        d2["padding"],
                        d2["compression"],
                        d2["crush_stress"],
                        d2["protection_area"],
                    )

                    energy_margin = capacity / max(
                        impact["impact_energy"], 1e-9
                    )

                    geometry_ok = (
                        d2["egg_length"] + 2 * d2["padding"]
                        <= d2["container_length"]
                        and d2["egg_width"] + 2 * d2["padding"]
                        <= d2["container_width"]
                        and d2["egg_height"] + 2 * d2["padding"]
                        <= d2["container_height"]
                    )

                    score = (
                        impact["heuristic_index"]
                        + max(0, 1 - min(energy_margin, 1)) * 60
                        + (50 if not geometry_ok else 0)
                    )

                                                                          
                                                                          
                                                                   
                                                                    
                    if not sim["landed"]:
                        score += 1000

                    rows.append(
                        (
                            score,
                            padding_mm,
                            structure_g,
                            sim["impact_velocity"],
                            impact["g_load"],
                            energy_margin,
                            geometry_ok,
                            sim["landed"],
                        )
                    )

            rows.sort(key=lambda x: x[0])

            self.optimization_text.delete("1.0", "end")
            self.optimization_text.insert(
                "end",
                "🔎 MELHORES CONFIGURAÇÕES ENCONTRADAS\n\n"
                "Score menor = combinação mais favorável segundo "
                "o modelo heurístico.\n\n"
                "Pad.     Massa     Impacto    G-load    Margem     Geometria\n"
                "─" * 72 + "\n",
            )

            for row in rows[:25]:
                score, padding, mass, velocity, g_load, margin, geometry, landed = row
                geometry_label = "OK" if geometry else "NÃO CABE"
                if not landed:
                    geometry_label += " (NÃO CONVERGIU)"
                self.optimization_text.insert(
                    "end",
                    f"{padding:5.0f} mm  "
                    f"{mass:6.1f} g  "
                    f"{velocity:7.2f} m/s  "
                    f"{g_load:7.1f} g  "
                    f"{margin:7.2f} x  "
                    f"{geometry_label}\n",
                )

            self.notebook.select(self.tab_save)

        except Exception as exc:
            messagebox.showerror("Erro na otimização", str(exc))

                                                              
            
                                                              

    def plot_placeholder(self):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.set_facecolor(PANEL)
        ax.text(
            0.5,
            0.5,
            "Execute uma simulação para gerar os gráficos.",
            ha="center",
            va="center",
            color=MUTED,
            transform=ax.transAxes,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        self.figure.patch.set_facecolor(BG)
        self.canvas.draw()

    def plot_results(self):
        if not self.result:
            return

        sim = self.result["simulation"]
        t = sim["time"]

        self.figure.clear()
        self.figure.patch.set_facecolor(BG)

        axes = self.figure.subplots(2, 2)
        titles = [
            "Velocidade",
            "Aceleração",
            "Altura",
            "Energia cinética",
        ]

        data = [
            (sim["velocity"], "m/s"),
            (sim["acceleration"], "m/s²"),
            (sim["height"], "m"),
            (sim["energy"], "J"),
        ]

        for ax, title, (values, ylabel) in zip(
            axes.flat, titles, data
        ):
            ax.set_facecolor(PANEL)
            ax.plot(t, values, linewidth=1.8)
            ax.set_title(title, color=TEXT)
            ax.set_xlabel("Tempo (s)", color=MUTED)
            ax.set_ylabel(ylabel, color=MUTED)
            ax.tick_params(colors=MUTED)
            for spine in ax.spines.values():
                spine.set_color("#334155")
            ax.grid(alpha=0.15)

        if sim["opening_time"] is not None:
            for ax in axes.flat:
                ax.axvline(
                    sim["opening_time"],
                    linestyle="--",
                    linewidth=1,
                )
                ylim = ax.get_ylim()
                ax.text(
                    sim["opening_time"],
                    ylim[1],
                    " abertura",
                    color=PURPLE,
                    va="top",
                )

        self.figure.tight_layout(pad=2)
        self.canvas.draw()

                                                              
                      
                                                              

    def project_dict(self):
        return {
            "version": "3.1",
            "mode": self.mode.get(),
            "variables": {
                key: value.get() for key, value in self.vars.items()
            },
        }

    def save_project(self):
        path = filedialog.asksaveasfilename(
            title="Salvar projeto",
            defaultextension=".json",
            filetypes=[("Projeto Egg Drop", "*.json"), ("Todos", "*.*")],
        )
        if not path:
            return

        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(
                    self.project_dict(),
                    f,
                    ensure_ascii=False,
                    indent=2,
                )
            self.last_project_path = path
            messagebox.showinfo(
                "Projeto salvo",
                f"Projeto salvo em:\n{path}",
            )
        except OSError as exc:
            messagebox.showerror("Erro ao salvar", str(exc))

    def load_project(self):
        path = filedialog.askopenfilename(
            title="Carregar projeto",
            filetypes=[("Projeto Egg Drop", "*.json"), ("Todos", "*.*")],
        )
        if not path:
            return

        try:
            with open(path, "r", encoding="utf-8") as f:
                project = json.load(f)

            variables = project.get("variables", {})
            for key, value in variables.items():
                if key in self.vars:
                    self.vars[key].set(str(value))

            if project.get("mode") in ("Básico", "Avançado"):
                self.mode.set(project["mode"])
            self._refresh_mode()

            self.last_project_path = path
            messagebox.showinfo(
                "Projeto carregado",
                "Os parâmetros do projeto foram carregados.",
            )

        except (OSError, json.JSONDecodeError) as exc:
            messagebox.showerror("Erro ao carregar", str(exc))

    def export_csv(self):
        if not self.result:
            messagebox.showwarning(
                "Nada para exportar",
                "Execute uma simulação antes de exportar.",
            )
            return

        path = filedialog.asksaveasfilename(
            title="Exportar resultados",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("Todos", "*.*")],
        )
        if not path:
            return

        d = self.result["data"]
        sim = self.result["simulation"]
        impact = self.result["impact"]

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
            ("capacidade_absorcao", self.result["absorption_capacity"], "J"),
            ("margem_energetica", self.result["energy_margin"], "x"),
            ("indice_heuristico", impact["heuristic_index"], "/100"),
            ("classificacao", impact["risk_label"], ""),
        ]

        try:
            with open(path, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f, delimiter=";")
                writer.writerows(rows)

            messagebox.showinfo(
                "CSV exportado",
                f"Resultados exportados para:\n{path}",
            )
        except OSError as exc:
            messagebox.showerror("Erro ao exportar", str(exc))

                                                              
                
                                                              

    def go_simulation(self):
        self.notebook.select(self.tab_sim)

    def go_parachute(self):
        self.notebook.select(self.tab_para)

    def go_protection(self):
        self.notebook.select(self.tab_protection)

    def go_graphs(self):
        self.notebook.select(self.tab_graphs)


def main():
    root = tk.Tk()
    app = EggDropLab(root)
    root.mainloop()


if __name__ == "__main__":
    main()
