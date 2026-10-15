#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Simplified theoretical level schemes:
plot lowest-energy calculated levels for each J^pi,
plus configurable extra states and strong branch transition arrows.
Target states of strong transitions are automatically highlighted,
and missing states involved in transitions are automatically included.

For:
    202Ir, 204Au, and 206Tl

Input:
    summary_jj56_Tl_Au_Ir.xlsx

Output:
    level_scheme_first_state/
        nucl_Ir202_first_level_each_Jpi.png / .pdf
        nucl_Au204_first_level_each_Jpi.png / .pdf
        nucl_Tl206_first_level_each_Jpi.png / .pdf
"""

from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt


# ==========================================================
# Settings
# ==========================================================

INPUT_XLSX = Path("summary_jj56_Tl_Au_Ir.xlsx")
OUTPUT_DIR = Path("level_scheme_first_state")

E_MAX = 3.0       # MeV
FIG_DPI = 300
SAVE_PDF = True

# Highlight target state
HIGHLIGHT_JPI = "12-"

# ==========================================================
# Manual half-life input (By Nucleus and Interaction Model)
# ==========================================================
HALF_LIFE_12MINUS = {
    "Ir202": {
        "jj56_ste_a0045": "T1/2 = 0.08147 s",
        "khhe_f_a01":      "T1/2 = 0.0233 s",
        "khhe_khhe_ste11": "T1/2 = 1.693 s",
    },
    "Au204": {
        "jj56_ste_a0045": "T1/2 = 0.05706 s\nH-like T1/2 = 1.259 s",
        "khhe_f_a01":      "T1/2 = 0.002146 s\nH-like T1/2 = 0.06386 s",
        "khhe_khhe_ste11": "T1/2 = 0.1223 s\nH-like T1/2 = 2.487 s",
    },
    "Tl206": {
        "jj56_ste_a0045": "T1/2 = 98.24 s",
        "khhe_f_a01":      "T1/2 = 98.82 s",
        "khhe_khhe_ste11": "T1/2 = 97.85 s",
    },
}

# ==========================================================
# Extra states configuration (Placeholder for second/third levels)
# ==========================================================
EXTRA_STATES = [
    # 示例占位：按需取消注释或添加您需要额外显示的能级
    # {"Nuclei": "Au204", "Int": "jj56_ste_a0045", "Jpi_plot": "8+", "N_level": 2},
    {"Nuclei": "Tl206", "Int": "jj56_ste_a0045", "Jpi_plot": "8+", "N_level": 2},
    {"Nuclei": "Tl206", "Int": "khhe_f_a01", "Jpi_plot": "8+", "N_level": 2},
    {"Nuclei": "Tl206", "Int": "khhe_khhe_ste11", "Jpi_plot": "8+", "N_level": 2},
]

# ==========================================================
# Strongest Branch Transitions Configuration
# 定义每个核素在各个模型下需要画箭头的强分支跃迁
# 格式: (Nucleus, Interaction): [(Initial_Jpi, Initial_N, Final_Jpi, Final_N), ...]
# ==========================================================
TRANSITIONS = {
    ("Au204", "jj56_ste_a0045"): [
        ("12-", 1, "8+", 1),
    ],
    ("Au204", "khhe_f_a01"): [
        ("12-", 1, "9+", 1),
    ],
    ("Au204", "khhe_khhe_ste11"): [
        ("12-", 1, "8+", 1),
    ],
    ("Ir202", "jj56_ste_a0045"): [
        ("12-", 1, "9+", 1),
    ],
    ("Ir202", "khhe_f_a01"): [
        ("12-", 1, "9+", 2),
    ],
    ("Ir202", "khhe_khhe_ste11"): [
        ("12-", 1, "8+", 1),
    ], 
    ("Tl206", "jj56_ste_a0045"): [
        ("12-", 1, "8+", 2),
    ],
    ("Tl206", "khhe_f_a01"): [
        ("12-", 1, "8+", 2),
    ],
    ("Tl206", "khhe_khhe_ste11"): [
        ("12-", 1, "8+", 2),
    ],      
}


OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

INTERACTIONS = [
    ("jj56_ste_a0045", "jj56_a01_h11g9_s1_ste_a0045"),
    ("khhe_f_a01", "khhe_f_a01"),
    ("khhe_khhe_ste11", "khhe_ste11"),
]

NUCLEI = [
    ("Ir202", r"$^{202}$Ir"),
    ("Au204", r"$^{204}$Au"),
    ("Tl206", r"$^{206}$Tl"),
]


# ==========================================================
# Read level data
# ==========================================================

def normalize_jpi(value):
    if pd.isna(value):
        return "?"
    return str(value).strip().replace(" ", "").replace("−", "-")


def load_levels():
    if not INPUT_XLSX.exists():
        raise FileNotFoundError(f"Cannot find {INPUT_XLSX.resolve()}")

    df = pd.read_excel(INPUT_XLSX, sheet_name="sheet1")
    required = {"Nuclei", "Int", "Jpi", "Ex", "E", "N_level"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    df = df.copy()
    df["Nuclei"] = df["Nuclei"].astype(str).str.strip()
    df["Int"] = df["Int"].astype(str).str.strip()
    df["Jpi_plot"] = df["Jpi"].map(normalize_jpi)

    for col in ["Ex", "E", "N_level"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["Ex", "E", "N_level"])
    df = df.drop_duplicates(subset=["Nuclei", "Int", "Jpi_plot", "Ex", "E", "N_level"])
    return df


# ==========================================================
# Select lowest state + extra states + transition related states
# ==========================================================

def get_selected_levels(sub, nucleus_key, interaction_key):
    sub_sorted = sub.sort_values(["Ex", "N_level", "E"])
    first_df = sub_sorted.groupby("Jpi_plot", sort=False, as_index=False).first()
    
    selected_frames = [first_df]
    
    # 1. 自动收集 TRANSITIONS 中涉及的能级并补上
    trans_list = TRANSITIONS.get((nucleus_key, interaction_key), [])
    trans_matches = []
    for (init_jpi, init_n, final_jpi, final_n) in trans_list:
        for jpi, n in [(init_jpi, init_n), (final_jpi, final_n)]:
            matched = sub[
                (sub["Jpi_plot"] == jpi) &
                (sub["N_level"] == n)
            ]
            if not matched.empty:
                trans_matches.append(matched)
                
    # 2. 收集 EXTRA_STATES 配置的能级
    if EXTRA_STATES:
        for rule in EXTRA_STATES:
            if rule.get("Nuclei") == nucleus_key and rule.get("Int") in [None, interaction_key]:
                matched = sub[
                    (sub["Jpi_plot"] == rule.get("Jpi_plot")) &
                    (sub["N_level"] == rule.get("N_level"))
                ]
                if not matched.empty:
                    trans_matches.append(matched)

    if trans_matches:
        extra_df = pd.concat(trans_matches, ignore_index=True)
        selected_frames.append(extra_df)
            
    combined = pd.concat(selected_frames, ignore_index=True)
    return combined.drop_duplicates(subset=["Jpi_plot", "N_level", "Ex"]).sort_values("Ex")


# ==========================================================
# Plot panel
# ==========================================================

def plot_panel(ax, levels, nucleus_key, interaction_key, interaction_label):
    sub = levels[
        (levels["Nuclei"] == nucleus_key) &
        (levels["Int"] == interaction_key) &
        (levels["Ex"] >= 0)
    ].copy()

    if sub.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    plot_levels = get_selected_levels(sub, nucleus_key, interaction_key)
    plot_levels = plot_levels[plot_levels["Ex"] <= E_MAX].sort_values("Ex")

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.05, E_MAX + 0.15)
    ax.set_xticks([])
    ax.grid(axis="y", linestyle=":", linewidth=0.4, alpha=0.35)
    ax.tick_params(axis="y", labelsize=8, length=0)

    for side in ["top", "right", "bottom", "left"]:
        ax.spines[side].set_visible(False)

    ax.set_title(interaction_label, fontsize=8.5, pad=10)

    # 找出当前模型下所有作为强跃迁终态的 (jpi, n_level) 集合
    trans_list = TRANSITIONS.get((nucleus_key, interaction_key), [])
    target_states = {(final_jpi, final_n) for (_, _, final_jpi, final_n) in trans_list}

    label_positions = []
    level_coords = {}

    for _, row in plot_levels.iterrows():
        y = float(row["Ex"])
        jpi = row["Jpi_plot"]
        n_lvl = int(row["N_level"])

        # 判断是否为高亮态（12- 或者是强跃迁的终态）
        is_highlight = (jpi.replace("−", "-") == HIGHLIGHT_JPI) or ((jpi, n_lvl) in target_states)

        color = "#1769aa" if is_highlight else "black"
        linewidth = 2.5 if is_highlight else 1.0

        ax.hlines(y, 0.30, 0.70, color=color, linewidth=linewidth, zorder=3)

        # 记录能级坐标
        level_coords[(jpi, n_lvl)] = y

        y_label = y
        for old in label_positions:
            if abs(y_label - old) < 0.06:
                y_label += 0.035
        label_positions.append(y_label)

        # J^pi 标签（如果是多重态，带上序号以区分，例如 8-[2] 或直接显示标签）
        label_text = f"{jpi}" if n_lvl == 1 else f"{jpi}({n_lvl})"

        ax.text(
            0.25, y_label, label_text,
            ha="right", va="center", fontsize=8,
            color=color, fontweight="bold" if is_highlight else "normal",
            clip_on=False
        )

        # 激发能数值
        ax.text(
            0.75, y_label, f"{y:.3f}",
            ha="left", va="center", fontsize=7.3,
            color=color, fontweight="bold" if is_highlight else "normal",
            clip_on=False
        )

        # 半衰期标注（仅对 12- 目标态显示）
        if jpi.replace("−", "-") == HIGHLIGHT_JPI and n_lvl == 1:
            model_halflives = HALF_LIFE_12MINUS.get(nucleus_key, {})
            half = model_halflives.get(interaction_key, "T1/2 = --- s")
            
            ax.text(
                0.86, y_label, f"({half})",
                ha="left", va="center", fontsize=6.8,
                color="#1769aa", fontweight="bold",
                clip_on=False
            )

    # ==========================================================
    # Draw Transition Arrows (Strongest branch transitions)
    # ==========================================================
    for (init_jpi, init_n, final_jpi, final_n) in trans_list:
        if (init_jpi, init_n) in level_coords and (final_jpi, final_n) in level_coords:
            y_init = level_coords[(init_jpi, init_n)]
            y_final = level_coords[(final_jpi, final_n)]
            
            ax.annotate(
                "",
                xy=(0.50, y_final), xycoords="data",
                xytext=(0.50, y_init), textcoords="data",
                arrowprops=dict(
                    arrowstyle="-|>",
                    color="#d9534f",      # 红色箭头
                    lw=1.5,
                    mutation_scale=8,
                    zorder=4
                )
            )


# ==========================================================
# Whole figure
# ==========================================================

def make_figure(levels, nucleus_key, nucleus_label):
    fig, axes = plt.subplots(1, 3, figsize=(14, 9), sharey=True, gridspec_kw={"wspace": 0.45})

    for i, (int_key, int_label) in enumerate(INTERACTIONS):
        plot_panel(axes[i], levels, nucleus_key, int_key, int_label)

    axes[0].set_ylabel(r"Excitation energy $E_x$ (MeV)", fontsize=10)
    fig.suptitle(f"Simplified theoretical level scheme: {nucleus_label}", fontsize=14, y=0.97)

    fig.text(
        0.5, 0.025,
        r"Lowest state for each $J^\pi$ shown + selected states. "
        r"Blue: first $12^-$ & transition targets. Red arrows: strongest branch transitions.",
        ha="center", fontsize=8.5
    )

    fig.subplots_adjust(left=0.09, right=0.95, top=0.91, bottom=0.08)

    safe = re.sub(r"[^A-Za-z0-9_-]+", "_", nucleus_key).strip("_")
    stem = OUTPUT_DIR / f"nucl_{safe}_first_level_each_Jpi"

    fig.savefig(stem.with_suffix(".png"), dpi=FIG_DPI, bbox_inches="tight")
    if SAVE_PDF:
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")

    plt.close(fig)
    print(f"Saved {stem.with_suffix('.png')} and {stem.with_suffix('.pdf')}")


# ==========================================================
# Main
# ==========================================================

def main():
    levels = load_levels()
    print(f"Loaded {len(levels)} states")

    for nucleus_key, nucleus_label in NUCLEI:
        make_figure(levels, nucleus_key, nucleus_label)

    print("\nOutput:", OUTPUT_DIR.resolve())


if __name__ == "__main__":
    main()