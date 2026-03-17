#!/usr/bin/env python3
"""
Wheel of Life - Visualisation of all species/genera across multiple markers

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle
from matplotlib.lines import Line2D
import pandas as pd
from pathlib import Path as FilePath
import warnings
warnings.filterwarnings('ignore')

SUPERCLADE_COLORS = {
    "Animals":        "#4CAF50",
    "Archaeplastida": "#E91E63",
    "Stramenopiles":  "#FF9800",
    "Alveolata":      "#2196F3",
    "Other":          "#9E9E9E",
}

GROUP_SUPERCLADE = {
    "Actinopteri":"Animals","Actinopterygii":"Animals",
    "Elasmobranchii":"Animals","Chondrichthyes":"Animals",
    "Amphibia":"Animals","Reptilia":"Animals","Aves":"Animals","Mammalia":"Animals",
    "Echinodermata":"Animals","Hemichordata":"Animals","Tunicata":"Animals",
    "Arthropoda":"Animals","Nematoda":"Animals","Tardigrada":"Animals",
    "Annelida":"Animals","Mollusca":"Animals","Bryozoa":"Animals","Brachiopoda":"Animals",
    "Platyhelminthes":"Animals","Nemertea":"Animals","Rotifera":"Animals","Sipuncula":"Animals",
    "Cnidaria":"Animals","Ctenophora":"Animals","Porifera":"Animals","Placozoa":"Animals",
    "Rhodophyta":"Archaeplastida","Chlorophyta":"Archaeplastida",
    "Streptophyta":"Archaeplastida","Plantae":"Archaeplastida","Charophyta":"Archaeplastida",
    "Bacillariophyta":"Stramenopiles","Ochrophyta":"Stramenopiles","Phaeophyta":"Stramenopiles",
    "Chrysophyta":"Stramenopiles","Xanthophyta":"Stramenopiles","Oomycota":"Stramenopiles",
    "Haptophyta":"Stramenopiles",
    "Dinoflagellata":"Alveolata","Dinophyceae":"Alveolata",
    "Ciliophora":"Alveolata","Apicomplexa":"Alveolata","Myzozoa":"Alveolata",
    "Cryptophyta":"Other","Foraminifera":"Other","Radiolaria":"Other","Cercozoa":"Other",
    "Euglenozoa":"Other","Heterolobosea":"Other","Metamonada":"Other",
    "Chytridiomycota":"Other","Basidiomycota":"Other","Ascomycota":"Other",
    "Choanoflagellatea":"Other","Pseudomonadota":"Other",
    "Other Major Group":"Other","Others":"Other","Other":"Other","Unknown":"Other",
}

# Phylogenetic order — groups present in the data but not listed fall back to size order
PHYLO_ORDER = [
    "Actinopteri","Actinopterygii","Elasmobranchii","Chondrichthyes",
    "Amphibia","Reptilia","Aves","Mammalia",
    "Echinodermata","Hemichordata","Tunicata",
    "Arthropoda","Nematoda","Tardigrada",
    "Annelida","Mollusca","Bryozoa","Brachiopoda",
    "Platyhelminthes","Nemertea","Rotifera","Sipuncula",
    "Cnidaria","Ctenophora","Porifera","Placozoa",
    "Chytridiomycota","Basidiomycota","Ascomycota","Choanoflagellatea",
    "Rhodophyta","Chlorophyta","Streptophyta","Plantae","Charophyta",
    "Bacillariophyta","Ochrophyta","Phaeophyta","Chrysophyta","Xanthophyta","Oomycota",
    "Dinoflagellata","Dinophyceae","Ciliophora","Apicomplexa","Myzozoa",
    "Foraminifera","Radiolaria","Cercozoa",
    "Haptophyta","Cryptophyta",
    "Euglenozoa","Heterolobosea","Metamonada",
    "Other Major Group",
]

# Shade palettes per superclade — assigned in PHYLO_ORDER to groups present in a wheel
_SHADES = {
    "Animals": [
        "#C8E6C9","#A5D6A7","#81C784","#66BB6A","#4CAF50","#43A047","#388E3C",
        "#2E7D32","#26A69A","#00897B","#00796B","#00695C","#004D40","#1B5E20",
        "#33691E","#558B2F","#689F38","#7CB342","#9CCC65","#AED581","#C5E1A5",
        "#DCEDC8","#F1F8E9","#E8F5E9",
    ],
    "Archaeplastida": [
        "#FCE4EC","#F8BBD9","#F48FB1","#F06292","#E91E63","#D81B60","#C2185B","#880E4F",
    ],
    "Stramenopiles": [
        "#FFF3E0","#FFE0B2","#FFCC80","#FFB74D","#FFA726","#FF9800","#FB8C00","#F57C00","#E65100",
    ],
    "Alveolata": [
        "#E3F2FD","#BBDEFB","#90CAF9","#64B5F6","#42A5F5","#2196F3","#1E88E5","#1565C0","#0D47A1",
    ],
    "Other": [
        "#FAFAFA","#F5F5F5","#EEEEEE","#E0E0E0","#BDBDBD","#9E9E9E","#757575","#616161","#424242",
    ],
}

_GROUP_COLOR_CACHE: dict = {}

# OTHERS left-to-right colour order: grey -> blue -> orange -> pink -> green
# Achieved by sorting OTHERS groups ascending by superclade rank (Other=0, Animals=4),
# so grey sorts first (leftmost) and green sorts last (rightmost).
# Colour shades: unified pass over ALL groups so OTHERS get the same shades as main.

_OTHERS_SC_RANK = {"Other": 0, "Alveolata": 1, "Stramenopiles": 2, "Archaeplastida": 3, "Animals": 4}

def _others_key(g, phylo_rank, n_phylo):
    """Sort key for OTHERS. Low = leftmost, high = rightmost.
    Primary:   superclade rank (Other=0 -> left/grey, Animals=4 -> right/green).
    Secondary: REVERSED phylo rank within each block so higher-phylo groups
               (which get darker shades from the unified counter) land on the
               LEFT, and lower-phylo (lighter shade) groups land on the RIGHT.
               Result: each superclade block goes dark-left -> light-right."""
    sc = GROUP_SUPERCLADE.get(g, "Other")
    BIG = n_phylo + 2
    sc_rank = _OTHERS_SC_RANK.get(sc, 0)
    within = BIG - 1 - phylo_rank.get(g, n_phylo)  # reversed: higher phylo -> smaller -> leftmost
    return sc_rank * BIG + within

def _build_color_cache(main_groups: list, others_groups: list = None):
    """Single unified PHYLO_ORDER pass over ALL groups (main + OTHERS).
    Every group gets the same shade it would get as a standalone main group."""
    _GROUP_COLOR_CACHE.clear()
    all_groups = list(main_groups) + (list(others_groups) if others_groups else [])
    ordered = [g for g in PHYLO_ORDER if g in all_groups] + \
              [g for g in all_groups if g not in PHYLO_ORDER]
    counters = {sc: 0 for sc in _SHADES}
    for g in ordered:
        sc = GROUP_SUPERCLADE.get(g, "Other")
        _GROUP_COLOR_CACHE[g] = _SHADES[sc][counters[sc] % len(_SHADES[sc])]
        counters[sc] += 1

def get_color(group):
    return _GROUP_COLOR_CACHE.get(group, "#9E9E9E")


# Marker symbols — 25 entries (5 shapes × 5 colours)
_MS_SHAPES  = ['o', 's', '^', 'D', 'v']
_MS_COLORS  = ['black', '#e41a1c', '#377eb8', '#4daf4a', '#984ea3']
MARKER_SYMBOLS = {}
_idx = 0
for _shape in _MS_SHAPES:
    for _color in _MS_COLORS:
        MARKER_SYMBOLS[_idx] = (_shape, _color, _color)
        _idx += 1


def find_contiguous_runs(positions):
    """
    Given a list of integer positions, find contiguous runs.
    Returns list of (start_pos, end_pos) tuples for each run.
    
    Example: [0,1,2,5,6,8] -> [(0,2), (5,6), (8,8)]
    """
    if not positions:
        return []
    
    positions = sorted(positions)
    runs = []
    start = positions[0]
    prev = positions[0]
    
    for pos in positions[1:]:
        if pos == prev + 1:
            prev = pos
        else:
            runs.append((start, prev))
            start = pos
            prev = pos
    
    runs.append((start, prev))
    return runs


def combine_marker_data(marker_dict):
    """
    Combine data from multiple markers.
    
    Parameters:
        marker_dict: dict of {marker_name: dataframe} e.g. {'12S': df_12s, '18S': df_18s}
    
    Returns:
        Combined dataframe with 'markers' column listing which markers detected each species per site
    """
    all_data = []
    marker_names = list(marker_dict.keys())
    
    for marker_name, df in marker_dict.items():
        df = df.copy()
        df['_marker'] = marker_name
        all_data.append(df)
    
    combined = pd.concat(all_data, ignore_index=True)

    # Normalise Class synonyms before groupby so Copepoda/Hexanauplia
    # cross-marker conflicts are resolved before first_valid() runs.
    _SYN = {'Copepoda': 'Hexanauplia', 'Maxillopoda': 'Hexanauplia'}
    if 'Class' in combined.columns:
        combined['Class'] = combined['Class'].replace(_SYN)

    # Columns to preserve (taxonomy + metadata + OTU counts)
    tax_cols = [c for c in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'display_group'] 
                if c in combined.columns]
    meta_cols = [c for c in ['Collector', 'Location', 'Date', 'n_otus', 'n_otu_groups'] 
                 if c in combined.columns]
    preserve_cols = tax_cols + meta_cols
    
    # Custom aggregation: take first non-null value
    def first_valid(series):
        valid = series.dropna()
        if len(valid) > 0:
            return valid.iloc[0]
        return None
    
    # For OTU counts, we need to SUM across markers (each marker contributes OTUs)
    # First get the site-level OTU counts from each marker
    if 'n_otus' in combined.columns:
        site_otu_counts = combined.groupby('Site').agg({
            'n_otus': 'sum',
            'n_otu_groups': 'max'  # max because groups overlap
        }).reset_index()
    else:
        site_otu_counts = None
    
    # For each Site+label combination, get unique markers and keep first valid taxonomy/metadata
    agg_dict = {'_marker': lambda x: ','.join(sorted(set(x)))}
    for col in preserve_cols:
        if col not in ['n_otus', 'n_otu_groups']:
            agg_dict[col] = first_valid
    
    grouped = combined.groupby(['Site', 'label']).agg(agg_dict).reset_index()
    grouped = grouped.rename(columns={'_marker': 'markers'})
    
    # Add back the summed OTU counts
    if site_otu_counts is not None:
        grouped = grouped.drop(columns=['n_otus', 'n_otu_groups'], errors='ignore')
        grouped = grouped.merge(site_otu_counts, on='Site', how='left')
    
    # Also need to preserve site-level metadata (Collector, Location, Date)
    # Get first valid per site from original combined data
    for col in ['Collector', 'Location', 'Date']:
        if col in combined.columns:
            site_meta = combined.groupby('Site')[col].agg(first_valid).reset_index()
            site_meta = site_meta.rename(columns={col: f'_site_{col}'})
            grouped = grouped.merge(site_meta, on='Site', how='left')
            # Use site-level if species-level is missing
            if col in grouped.columns:
                grouped[col] = grouped[col].fillna(grouped[f'_site_{col}'])
            else:
                grouped[col] = grouped[f'_site_{col}']
            grouped = grouped.drop(columns=[f'_site_{col}'])
    
    print(f"Combined {len(marker_dict)} markers:")
    for name in marker_names:
        n = len(marker_dict[name])
        cols = list(marker_dict[name].columns)
        print(f"  {name}: {n} records, columns: {cols}")
    print(f"  Total unique species: {len(grouped)} records")
    print(f"  Final columns: {list(grouped.columns)}")
    
    return grouped, marker_names


def load_and_combine_markers(file_dict):
    """
    Load CSV files and combine them.
    
    Parameters:
        file_dict: dict of {marker_name: filepath} e.g. {'12S': 'wheel_data_12S.csv', '18S': 'wheel_data_18S.csv'}
    
    Returns:
        Combined dataframe with marker info
    """
    marker_dict = {}
    for name, filepath in file_dict.items():
        # Use latin-1 encoding to handle special characters
        marker_dict[name] = pd.read_csv(filepath, encoding='latin-1')
    
    return combine_marker_data(marker_dict)


def get_taxonomy_columns(df):
    """Get available taxonomy columns in hierarchical order."""
    tax_order = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    available = [col for col in tax_order if col in df.columns]
    return available


def draw_hierarchical_tree(ax, df, angles_rad, tree_color, root_r, tip_r):
    """
    Bifurcating cladogram: root → Class → Order → Family → Genus → tip.

    CRITICAL DESIGN CHOICES that prevent visual mis-grouping:

    1. ALL midpoint/span maths use POSITION-INDEX space (0…n-1), a clean
       linear sequence. Radian conversion only happens at the final drawing
       step. This eliminates angular wraparound artefacts near the 90° seam.

    2. SIX distinct radius levels — each taxonomic level occupies its own
       ring, so arcs from different levels NEVER overlap:
           root → R_class → R_order → R_family → R_genus → R_species → tip
       R_genus  = spanning arc connecting different genera within one family.
       R_species = spanning arc connecting different species within one genus.
       Without this separation the family connector arc sat at the same radius
       as the genus arcs, making adjacent genera appear grouped together.

    3. At every level: ONE spanning arc from leftmost to rightmost child mid,
       plus a spoke from the parent's mid-position to the mid of that span.
    """
    n_species = len(df)
    if n_species == 0:
        return

    R_class   = root_r + 0.025
    R_order   = root_r + 0.065
    R_family  = root_r + 0.115
    R_genus   = tip_r  - 0.055  # family spanning arc
    R_species = tip_r  - 0.025  # genus spanning arc

    group_col = 'original_display_group' if 'original_display_group' in df.columns else 'display_group'

    _SYN = {'Copepoda': 'Hexanauplia', 'Maxillopoda': 'Hexanauplia'}
    df = df.copy()
    for col in ['Class', 'Order', 'Family', 'Genus']:
        if col in df.columns:
            df[col] = df[col].fillna('Unknown')
    if 'Class' in df.columns:
        df['Class'] = df['Class'].replace(_SYN)

    # --- helpers (position space only) -----------------------------------

    def _pos(idx_list):
        return sorted([df.index.get_loc(i) for i in idx_list])

    def _pmid(positions):
        return (min(positions) + max(positions)) / 2.0

    def _p2r(pos):
        return np.radians(90.0 - (pos / n_species) * 360.0)

    def _arc(r, p_lo, p_hi, lw):
        if abs(p_hi - p_lo) < 0.01:
            return
        a1, a2 = _p2r(p_hi), _p2r(p_lo)   # higher pos → lower angle
        pts = np.linspace(a1, a2, max(3, int(50 * abs(a2 - a1) / np.pi)))
        ax.plot(r * np.cos(pts), r * np.sin(pts),
                color=tree_color, linewidth=lw, zorder=4)

    def _spoke(r1, r2, pos, lw):
        a = _p2r(pos)
        ax.plot([r1 * np.cos(a), r2 * np.cos(a)],
                [r1 * np.sin(a), r2 * np.sin(a)],
                color=tree_color, linewidth=lw, zorder=4, solid_capstyle='round')

    def _arc_to_children(r, parent_pos, child_mids, lw):
        """Arc at radius r spanning from parent_pos through all child_mids.
        parent_pos is guaranteed ON the arc; each child_mid is also on it.
        The incoming spoke always touches at parent_pos; outgoing spokes
        each depart at their child_mid — all on the same arc, no gaps."""
        all_pts = [parent_pos] + list(child_mids)
        lo, hi = min(all_pts), max(all_pts)
        if abs(hi - lo) < 0.005:
            return
        a1, a2 = _p2r(hi), _p2r(lo)
        pts = np.linspace(a1, a2, max(3, int(60 * abs(a2 - a1) / np.pi)))
        ax.plot(r * np.cos(pts), r * np.sin(pts),
                color=tree_color, linewidth=lw, zorder=4, solid_capstyle='round')

    # --- tree ------------------------------------------------------------

    for dg in df[group_col].unique():
        dg_idx = df[df[group_col] == dg].index.tolist()
        if not dg_idx:
            continue
        if len(dg_idx) == 1:
            _spoke(root_r, tip_r, _pos(dg_idx)[0], 0.5)
            continue

        # Each node's position = midpoint of ALL its descendent species (not just child mids).
        # This ensures the parent spoke always arrives exactly where the child arc starts.

        # ── Class ────────────────────────────────────────────────────────
        has_class = 'Class' in df.columns
        cls_info = {}   # cls → (mid_of_all_cls_species, [idx])
        for cls in (df.loc[dg_idx, 'Class'].unique() if has_class else ['Unknown']):
            c_idx = df[(df[group_col] == dg) & (df['Class'] == cls)].index.tolist() \
                    if has_class else dg_idx
            if c_idx:
                cls_info[cls] = (_pmid(_pos(c_idx)), c_idx)

        dg_mid  = _pmid(_pos(dg_idx))          # group root position
        cls_mids = [v[0] for v in cls_info.values()]
        _spoke(root_r, R_class, dg_mid, 0.8)
        _arc_to_children(R_class, dg_mid, cls_mids, 0.7)

        for cls, (cls_mid, c_idx) in cls_info.items():
            # cls_mid = midpoint of ALL species in this class → spoke departs here
            if len(c_idx) == 1:
                _spoke(R_class, tip_r, _pos(c_idx)[0], 0.5)
                continue

            # ── Order ────────────────────────────────────────────────────
            has_order = 'Order' in df.columns
            ord_info = {}
            for order in (df.loc[c_idx, 'Order'].unique() if has_order else ['Unknown']):
                o_idx = df[(df[group_col] == dg) & (df['Class'] == cls) &
                           (df['Order'] == order)].index.tolist() if has_order else c_idx
                if o_idx:
                    ord_info[order] = (_pmid(_pos(o_idx)), o_idx)

            ord_mids = [v[0] for v in ord_info.values()]
            _spoke(R_class, R_order, cls_mid, 0.6)      # depart from cls_mid ← key fix
            _arc_to_children(R_order, cls_mid, ord_mids, 0.6)

            for order, (ord_mid, o_idx) in ord_info.items():
                if len(o_idx) == 1:
                    _spoke(R_order, tip_r, _pos(o_idx)[0], 0.5)
                    continue

                # ── Family ───────────────────────────────────────────────
                has_family = 'Family' in df.columns
                fam_info = {}
                for fam in (df.loc[o_idx, 'Family'].unique() if has_family else ['Unknown']):
                    f_idx = df[(df[group_col] == dg) & (df['Class'] == cls) &
                               (df['Order'] == order) &
                               (df['Family'] == fam)].index.tolist() if has_family else o_idx
                    if f_idx:
                        fam_info[fam] = (_pmid(_pos(f_idx)), f_idx)

                fam_mids = [v[0] for v in fam_info.values()]
                _spoke(R_order, R_family, ord_mid, 0.5)     # depart from ord_mid ← key fix
                _arc_to_children(R_family, ord_mid, fam_mids, 0.5)

                for fam, (fam_mid, f_idx) in fam_info.items():
                    if len(f_idx) == 1:
                        _spoke(R_family, tip_r, _pos(f_idx)[0], 0.4)
                        continue

                    # ── Genus ─────────────────────────────────────────────
                    has_genus = 'Genus' in df.columns
                    genus_runs = []
                    for genus in (df.loc[f_idx, 'Genus'].unique() if has_genus else ['Unknown']):
                        g_idx = df[(df[group_col] == dg) & (df['Class'] == cls) &
                                   (df['Order'] == order) & (df['Family'] == fam) &
                                   (df['Genus'] == genus)].index.tolist() if has_genus else f_idx
                        if not g_idx:
                            continue
                        for rs, re in find_contiguous_runs(
                                sorted([df.index.get_loc(i) for i in g_idx])):
                            run = list(range(rs, re + 1))
                            genus_runs.append((_pmid(run), run))

                    if not genus_runs:
                        continue

                    gen_mids = [r[0] for r in genus_runs]
                    _spoke(R_family, R_genus, fam_mid, 0.4)     # depart from fam_mid ← key fix
                    _arc_to_children(R_genus, fam_mid, gen_mids, 0.4)

                    for (gen_mid, run_pos) in genus_runs:
                        _spoke(R_genus, R_species, gen_mid, 0.35)  # depart from gen_mid ← key fix
                        _arc_to_children(R_species, gen_mid, run_pos, 0.3)
                        for p in run_pos:
                            _spoke(R_species, tip_r, p, 0.25)


def create_wheel_of_life(df, title="", output_path=None, figsize=14, show_tree=True, 
                         marker_names=None, show_legend=True,
                         collector=None, location=None, date=None,
                         center_image=None, center_image_zoom=0.04,
                         include_genus=True, min_group_size=3):
    
    df = df.dropna(subset=['label', 'display_group']).drop_duplicates(subset=['label'])
    n_species = len(df)
    
    if n_species < 3:
        print(f"Need at least 3 species, got {n_species}")
        return None, None
    
    # Group small display_groups into "Other Major Group" but keep original for coloring
    df = df.copy()
    df['original_display_group'] = df['display_group']  # Save original for colors
    group_counts = df['display_group'].value_counts()
    small_groups = group_counts[group_counts < min_group_size].index.tolist()
    if small_groups:
        df.loc[df['display_group'].isin(small_groups), 'display_group'] = 'Other Major Group'

    # Build colour cache: main display groups get standard shade assignment;
    # OTHERS original groups get fresh assignment from lightest shade.
    main_display_groups = [g for g in df['display_group'].unique() if g != 'Other Major Group']
    others_original_groups = list(df.loc[df['display_group'] == 'Other Major Group',
                                         'original_display_group'].unique())
    _build_color_cache(main_display_groups, others_original_groups)
    # Debug: print OTHERS group order so you can verify right→left = green→grey
    if others_original_groups:
        _pr = {g: i for i, g in enumerate(PHYLO_ORDER)}
        _np = len(PHYLO_ORDER)
        _sorted_dbg = sorted(others_original_groups,
                             key=lambda g: _others_key(g, _pr, _np))
        print(f"OTHERS sort order (clockwise/right→counterclockwise/left = green→grey): {_sorted_dbg}")

    # Get taxonomy columns
    tax_cols = get_taxonomy_columns(df)

    # Order display groups phylogenetically
    group_counts = df['display_group'].value_counts()
    present_groups = group_counts.index.tolist()
    phylo_ranked = [g for g in PHYLO_ORDER if g in present_groups and g != 'Other Major Group']
    size_ranked  = [g for g in present_groups if g not in PHYLO_ORDER and g != 'Other Major Group']
    group_order  = phylo_ranked + size_ranked
    if 'Other Major Group' in present_groups:
        group_order.append('Other Major Group')

    # Normalise Class synonyms before sort so same-family species land contiguously
    _SYN = {'Copepoda': 'Hexanauplia', 'Maxillopoda': 'Hexanauplia'}
    if 'Class' in df.columns:
        df['Class'] = df['Class'].replace(_SYN)

    df['group_order'] = df['display_group'].map({g: i for i, g in enumerate(group_order)})

    # OTHERS sub-sort: Animals(sc=0) first → clockwise/right; Other(sc=4) last → left.
    # Use single integer key _others_key to avoid pandas tuple-column sorting issues.
    phylo_rank = {g: i for i, g in enumerate(PHYLO_ORDER)}
    n_phylo = len(PHYLO_ORDER)
    df['_others_sort'] = df['original_display_group'].map(
        lambda g: _others_key(g, phylo_rank, n_phylo))

    sort_cols = ['group_order', '_others_sort'] + \
                [col for col in ['Class', 'Order', 'Family', 'Genus'] if col in df.columns] + ['label']
    df = df.sort_values([c for c in sort_cols if c in df.columns]).reset_index(drop=True)
    df = df.drop(columns=['_others_sort'], errors='ignore')
    
    n_groups = len(group_order)
    
    # Assign angular positions
    positions = np.arange(n_species)
    angles_deg = 90 - (positions / n_species) * 360
    df['angle_deg'] = angles_deg
    df['angle_rad'] = np.radians(angles_deg)
    
    angles_rad = df['angle_rad'].tolist()
    
    # Group info for wedges and outer ring
    groups_info = []
    for group in group_order:
        group_df = df[df['display_group'] == group]
        half_slot = 180 / n_species
        start = group_df['angle_deg'].max() + half_slot
        end = group_df['angle_deg'].min() - half_slot
        mid = (start + end) / 2
        groups_info.append({
            'group': group,
            'start_deg': start,
            'end_deg': end,
            'mid_deg': mid,
            'mid_rad': np.radians(mid),
            'color': get_color(group),
            'n': len(group_df),
        })
    
    # Build sub-wedge info for coloring (preserves original phylum colors within "Other Phyla")
    subwedges_info = []
    for group in group_order:
        group_df = df[df['display_group'] == group]
        half_slot = 180 / n_species
        
        if group == 'Other Major Group':
            orig_groups_present = group_df['original_display_group'].unique().tolist()
            orig_sorted = sorted(
                orig_groups_present,
                key=lambda g: _others_key(g, phylo_rank, n_phylo)
            )
            for orig_group in orig_sorted:
                orig_df = group_df[group_df['original_display_group'] == orig_group]
                if orig_df.empty:
                    continue
                start = orig_df['angle_deg'].max() + half_slot
                end = orig_df['angle_deg'].min() - half_slot
                subwedges_info.append({
                    'start_deg': start,
                    'end_deg': end,
                    'color': get_color(orig_group),
                })
        else:
            start = group_df['angle_deg'].max() + half_slot
            end = group_df['angle_deg'].min() - half_slot
            subwedges_info.append({
                'start_deg': start,
                'end_deg': end,
                'color': get_color(group),
            })
    
    fig, ax = plt.subplots(figsize=(figsize, figsize))
    ax.set_aspect('equal')
    ax.axis('off')
    
    # RADIAL ZONES
    root_r = 0.08
    tip_r = 0.28
    wedge_inner = 0.30       
    wedge_outer = 0.76      
    outer_ring_inner = 0.79
    outer_ring_outer = 0.92
    
    # 1. Colored wedges (using subwedges to preserve original phylum colors)
    for info in subwedges_info:
        wedge = Wedge((0, 0), wedge_outer, theta1=info['end_deg'], theta2=info['start_deg'],
                      width=wedge_outer - wedge_inner, facecolor=info['color'],
                      edgecolor='white', linewidth=1.5)
        ax.add_patch(wedge)
    
    # 2. Outer dark ring
    for info in groups_info:
        ring = Wedge((0, 0), outer_ring_outer, theta1=info['end_deg'], theta2=info['start_deg'],
                     width=outer_ring_outer - outer_ring_inner, facecolor='#1a237e',
                     edgecolor='white', linewidth=0.5)
        ax.add_patch(ring)
    
    # 3. White center
    ax.add_patch(Circle((0, 0), root_r - 0.01, facecolor='white', edgecolor='none', zorder=5))
    
    # 3b. CENTER IMAGE (if provided)
    if center_image is not None:
        try:
            from matplotlib.offsetbox import OffsetImage, AnnotationBbox
            import matplotlib.image as mpimg
            
            img = mpimg.imread(center_image)
            
            # Use the zoom parameter - smaller number = smaller image
            imagebox = OffsetImage(img, zoom=center_image_zoom)
            ab = AnnotationBbox(imagebox, (0, 0), frameon=False, zorder=6)
            ax.add_artist(ab)
        except Exception as e:
            print(f"Warning: Could not load center image: {e}")
    
    # 4. PHYLOGENETIC TREE
    if show_tree:
        tree_color = '#505050'
        
        # Central root arc
        theta_full = np.linspace(0, 2*np.pi, 200)
        ax.plot(root_r * np.cos(theta_full), root_r * np.sin(theta_full),
                color=tree_color, linewidth=1.0, zorder=4)
        
        # Draw hierarchical tree
        draw_hierarchical_tree(ax, df, angles_rad, tree_color, root_r, tip_r)
    
    # 5. SPECIES LABELS with marker symbols
    fontsize = 7
    if n_species > 50: fontsize = 5
    elif n_species > 35: fontsize = 6
    elif n_species < 15: fontsize = 9
    
    text_offset = 0.02
    symbol_size = 4 if n_species > 30 else 5
    
    # Build marker index lookup
    marker_idx = {}
    if marker_names:
        for i, name in enumerate(marker_names):
            marker_idx[name] = i
    
    for _, row in df.iterrows():
        angle_rad = row['angle_rad']
        angle_deg = row['angle_deg']
        
        norm_angle = angle_deg
        while norm_angle > 180: norm_angle -= 360
        while norm_angle < -180: norm_angle += 360
        
        # Estimate text width in data coordinates
        text_width = len(row['label']) * fontsize * 0.0013
        
        if -90 <= norm_angle <= 90:
            # RHS: text reads outward, symbols after last letter (outer end)
            rotation = angle_deg
            ha = 'left'
            text_r = wedge_inner + text_offset
        else:
            # LHS: text flipped, symbols before first letter (also outer end)
            rotation = angle_deg + 180
            ha = 'right'
            text_r = wedge_inner + text_offset
        
        x = text_r * np.cos(angle_rad)
        y = text_r * np.sin(angle_rad)
        
        # Draw species label
        ax.text(x, y, row['label'],
                rotation=rotation,
                rotation_mode='anchor',
                ha=ha, va='center',
                fontsize=fontsize,
                fontstyle='italic',
                color='black',
                zorder=10)
        
        # Draw marker symbols at outer end of text (both sides)
        if 'markers' in row.index and marker_names and pd.notna(row['markers']):
            species_markers = row['markers'].split(',')
            
            # Symbols go at outer end of text (text_r + text_width + small gap)
            symbol_base_r = text_r + text_width + 0.012
            
            for i, marker in enumerate(species_markers):
                marker = marker.strip()
                if marker in marker_idx:
                    idx = marker_idx[marker] % len(MARKER_SYMBOLS)
                    symbol, facecolor, edgecolor = MARKER_SYMBOLS[idx]
                    
                    # Position symbols along radial line, spaced outward
                    sym_r = symbol_base_r + i * 0.012
                    sym_x = sym_r * np.cos(angle_rad)
                    sym_y = sym_r * np.sin(angle_rad)
                    
                    ax.plot(sym_x, sym_y, marker=symbol, markersize=symbol_size,
                           markerfacecolor=facecolor, markeredgecolor=edgecolor,
                           markeredgewidth=0.8, zorder=11)
    
    # 6. GROUP LABELS
    group_label_r = (outer_ring_inner + outer_ring_outer) / 2
    
    for info in groups_info:
        text = 'OTHERS' if info['group'] == 'Other Major Group' else info['group'].upper()
        mid_deg = info['mid_deg']
        mid_rad = info['mid_rad']
        
        x = group_label_r * np.cos(mid_rad)
        y = group_label_r * np.sin(mid_rad)
        
        norm_angle = mid_deg
        while norm_angle > 180: norm_angle -= 360
        while norm_angle < -180: norm_angle += 360
        
        if 0 <= norm_angle <= 180:
            rotation = norm_angle - 90
        else:
            rotation = norm_angle + 90
        
        arc_span = abs(info['start_deg'] - info['end_deg'])
        
        # More aggressive scaling for small arcs
        # Base font size on arc span directly
        if arc_span < 5:
            grp_fontsize = 3
        elif arc_span < 10:
            grp_fontsize = 4
        elif arc_span < 15:
            grp_fontsize = 5
        elif arc_span < 25:
            grp_fontsize = 6
        elif arc_span < 40:
            grp_fontsize = 7
        else:
            grp_fontsize = min(10, max(6, 8 * arc_span / 50))
        
        ax.text(x, y, text, rotation=rotation, ha='center', va='center',
                fontsize=grp_fontsize, fontweight='bold', color='white', zorder=10)
    
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    
    # Get OTU counts (total before species filtering) if available
    if 'n_otus' in df.columns and pd.notna(df['n_otus'].iloc[0]):
        n_otus = int(df['n_otus'].iloc[0])
        n_otu_groups = int(df['n_otu_groups'].iloc[0]) if 'n_otu_groups' in df.columns else n_groups
    else:
        n_otus = n_species
        n_otu_groups = n_groups
    
    # n_species and n_groups are the species-level counts (what's shown on wheel)
    
    if title:
        ax.text(0, 0.97, title, ha='center', va='bottom', fontsize=14, fontweight='bold')
        taxa_label = "genera/species" if include_genus else "species"
        subtitle = f"{n_otus} OTUs in {n_otu_groups} groups identified to {n_species} {taxa_label} in {n_groups} groups"
        ax.text(0, 0.92, subtitle, ha='center', va='bottom', fontsize=9, color='gray')
    
    # 7. MARKER LEGEND (bottom right)
    if marker_names and show_legend and len(marker_names) > 0:
        legend_elements = []
        for i, name in enumerate(marker_names):
            idx = i % len(MARKER_SYMBOLS)
            symbol, facecolor, edgecolor = MARKER_SYMBOLS[idx]
            legend_elements.append(
                Line2D([0], [0], marker=symbol, color='w', 
                       markerfacecolor=facecolor, markeredgecolor=edgecolor,
                       markersize=8, markeredgewidth=1.0, label=name)
            )
        
        legend = ax.legend(handles=legend_elements, loc='lower right',
                          bbox_to_anchor=(0.98, 0.02),
                          fontsize=8, frameon=True, fancybox=True,
                          framealpha=0.9, edgecolor='gray',
                          title='Markers', title_fontsize=9)

    # 7b. SUPERCLADE LEGEND (top right)
    # Shows every superclade present in this wheel (including those that only appear
    # in the OTHERS wedge), each with its midpoint colour.
    # Uses an inset axes so it never clashes with the wheel or marker legend.
    sc_present = []
    sc_order = ["Animals", "Archaeplastida", "Stramenopiles", "Alveolata", "Other"]
    for sc in sc_order:
        # Check main groups first
        groups_in_sc = [g for g in group_order if GROUP_SUPERCLADE.get(g, "Other") == sc]
        if groups_in_sc:
            mid_g = groups_in_sc[len(groups_in_sc) // 2]
            sc_present.append((sc, _GROUP_COLOR_CACHE.get(mid_g, SUPERCLADE_COLORS[sc])))
        else:
            # Also check groups that were folded into OTHERS
            others_in_sc = [g for g in others_original_groups
                            if GROUP_SUPERCLADE.get(g, "Other") == sc]
            if others_in_sc:
                mid_g = others_in_sc[len(others_in_sc) // 2]
                sc_present.append((sc, _GROUP_COLOR_CACHE.get(mid_g, SUPERCLADE_COLORS[sc])))

    if sc_present:
        from matplotlib.patches import Patch
        sc_handles = [Patch(facecolor=col, edgecolor='white', linewidth=0.5, label=name)
                      for name, col in sc_present]
        sc_legend = ax.legend(handles=sc_handles, loc='upper right',
                              bbox_to_anchor=(0.98, 0.98),
                              fontsize=8, frameon=True, fancybox=True,
                              framealpha=0.9, edgecolor='gray',
                              title='Major Groups', title_fontsize=9,
                              handlelength=1.2, handleheight=1.2)
        # Re-add marker legend if it was drawn (adding a second legend replaces the first)
        if marker_names and show_legend and len(marker_names) > 0:
            ax.add_artist(sc_legend)
            ax.legend(handles=legend_elements, loc='lower right',
                      bbox_to_anchor=(0.98, 0.02),
                      fontsize=8, frameon=True, fancybox=True,
                      framealpha=0.9, edgecolor='gray',
                      title='Markers', title_fontsize=9)
    import unicodedata
    
    def clean_text(val):
        """Clean text for display, handling None/NA and special characters."""
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return "NA"
        text = str(val).strip()
        if text == "" or text.lower() == "nan" or text.lower() == "na":
            return "NA"
        
        # Normalize unicode and replace accented characters with ASCII equivalents
        # First try to normalize to decomposed form (NFD) then strip combining marks
        try:
            # Convert to NFD (decomposed) form
            normalized = unicodedata.normalize('NFD', text)
            # Keep only non-combining characters (removes accents)
            ascii_text = ''.join(c for c in normalized if unicodedata.category(c) != 'Mn')
            
            # Manual replacements for characters that don't decompose well
            replacements = {
                'ø': 'o', 'Ø': 'O',
                'æ': 'ae', 'Æ': 'AE',
                'ß': 'ss',
                'ð': 'd', 'Ð': 'D',
                'þ': 'th', 'Þ': 'Th',
                'œ': 'oe', 'Œ': 'OE',
                ''': "'", ''': "'",
                '"': '"', '"': '"',
                '–': '-', '—': '-',
                '…': '...',
            }
            for old, new in replacements.items():
                ascii_text = ascii_text.replace(old, new)
            
            # Final cleanup - encode to ASCII, replacing any remaining non-ASCII
            ascii_text = ascii_text.encode('ascii', 'replace').decode('ascii')
            return ascii_text
        except Exception:
            # Fallback: just replace non-ASCII with ?
            return text.encode('ascii', 'replace').decode('ascii')
    
    collector_clean = clean_text(collector)
    location_clean = clean_text(location)
    date_clean = clean_text(date)
    
    info_lines = [
        f"Collector: {collector_clean}",
        f"Location: {location_clean}",
        f"Date: {date_clean}"
    ]
    
    info_text = '\n'.join(info_lines)
    # Add text box in bottom left using data coordinates
    props = dict(boxstyle='round,pad=0.4', facecolor='white', 
                 edgecolor='gray', alpha=0.9)
    ax.text(-0.95, -0.88, info_text,
            fontsize=8, verticalalignment='top', horizontalalignment='left',
            bbox=props, zorder=15)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"Saved: {output_path}")
    
    return fig, ax


def generate_site_wheels(df, output_dir="Figures/wheels", prefix="", 
                         format="pdf", min_species=5, marker_names=None, 
                         include_genus=True, min_group_size=3, **kwargs):
    output_dir = FilePath(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Debug: show what columns we have
    print(f"\nColumns in data: {list(df.columns)}")
    for col in ['Collector', 'Location', 'Date']:
        if col in df.columns:
            non_null = df[col].notna().sum()
            print(f"  {col}: {non_null} non-null values")
    
    sites = df['Site'].unique()
    print(f"\n=== Generating {len(sites)} wheels ===")
    
    created = skipped = 0
    for site in sites:
        site_df = df[df['Site'] == site]
        n_sp = site_df['label'].nunique()
        
        if n_sp < min_species:
            skipped += 1
            continue
        
        print(f"{site} ({n_sp} sp)... ", end="", flush=True)
        safe_name = "".join(c if c.isalnum() or c in '_-' else '_' for c in site)
        file_prefix = f"{prefix}_" if prefix else ""
        filename = output_dir / f"{file_prefix}wheel_{safe_name}.{format}"
        
        # Extract metadata - find first non-null value in site
        collector = None
        location = None
        date = None
        
        if 'Collector' in site_df.columns:
            valid_vals = site_df['Collector'].dropna()
            valid_vals = valid_vals[valid_vals.astype(str).str.strip() != '']
            if len(valid_vals) > 0:
                collector = str(valid_vals.iloc[0]).strip()
        
        if 'Location' in site_df.columns:
            valid_vals = site_df['Location'].dropna()
            valid_vals = valid_vals[valid_vals.astype(str).str.strip() != '']
            if len(valid_vals) > 0:
                location = str(valid_vals.iloc[0]).strip()
        
        if 'Date' in site_df.columns:
            valid_vals = site_df['Date'].dropna()
            valid_vals = valid_vals[valid_vals.astype(str).str.strip() != '']
            if len(valid_vals) > 0:
                date = str(valid_vals.iloc[0]).strip()
        
        try:
            fig, ax = create_wheel_of_life(site_df, title=site, output_path=filename, 
                                           marker_names=marker_names,
                                           collector=collector, location=location, date=date,
                                           include_genus=include_genus, min_group_size=min_group_size,
                                           **kwargs)
            if fig:
                plt.close(fig)
                print("OK")
                created += 1
        except Exception as e:
            print(f"ERROR: {e}")
    
    print(f"\nCreated: {created}, Skipped: {skipped}")


if __name__ == "__main__":
    print("Creating demo wheel with hierarchical tree...")
    
    # -------------------------------------------------------------------------
    # Demo simulates how the real pipeline works:
    # - Each marker is a separate DataFrame (as extracted from phyloseq in R)
    # - combine_marker_data() assembles the 'markers' column dynamically from
    #   whatever keys you pass in — no hardcoded marker strings anywhere
    # - Change the dict keys below to match your own marker names
    # -------------------------------------------------------------------------
    
    # Shared taxonomy helper
    _base_tax = {
        'Dicentrarchus labrax': ('Actinopteri', 'Chordata',   'Actinopteri',  'Perciformes',   'Moronidae',    'Dicentrarchus'),
        'Sparus aurata':        ('Actinopteri', 'Chordata',   'Actinopteri',  'Spariformes',   'Sparidae',     'Sparus'),
        'Diplodus vulgaris':    ('Actinopteri', 'Chordata',   'Actinopteri',  'Spariformes',   'Sparidae',     'Diplodus'),
        'Sardina pilchardus':   ('Actinopteri', 'Chordata',   'Actinopteri',  'Clupeiformes',  'Clupeidae',    'Sardina'),
        'Chelon auratus':       ('Actinopteri', 'Chordata',   'Actinopteri',  'Mugiliformes',  'Mugilidae',    'Chelon'),
        'Mullus surmuletus':    ('Actinopteri', 'Chordata',   'Actinopteri',  'Mulliformes',   'Mullidae',     'Mullus'),
        'Patella caerulea':     ('Mollusca',    'Mollusca',   'Gastropoda',   'Patellida',     'Patellidae',   'Patella'),
        'Patella vulgata':      ('Mollusca',    'Mollusca',   'Gastropoda',   'Patellida',     'Patellidae',   'Patella'),
        'Octopus vulgaris':     ('Mollusca',    'Mollusca',   'Cephalopoda',  'Octopoda',      'Octopodidae',  'Octopus'),
        'Mytilus edulis':       ('Mollusca',    'Mollusca',   'Bivalvia',     'Mytilida',      'Mytilidae',    'Mytilus'),
        'Carcinus maenas':      ('Arthropoda',  'Arthropoda', 'Malacostraca', 'Decapoda',      'Portunidae',   'Carcinus'),
        'Balanus trigonus':     ('Arthropoda',  'Arthropoda', 'Hexanauplia',  'Balanomorpha',  'Balanidae',    'Balanus'),
        'Paracalanus parvus':   ('Arthropoda',  'Arthropoda', 'Hexanauplia',  'Calanoida',     'Paracalanidae','Paracalanus'),
        'Acartia clausii':      ('Arthropoda',  'Arthropoda', 'Hexanauplia',  'Calanoida',     'Acartiidae',   'Acartia'),
        'Paracentrotus lividus':('Echinodermata','Echinodermata','Echinoidea', 'Camarodonta',  'Parechinidae', 'Paracentrotus'),
        'Asterias rubens':      ('Echinodermata','Echinodermata','Asteroidea', 'Forcipulatida','Asteriidae',   'Asterias'),
        'Homo sapiens':         ('Mammalia',    'Chordata',   'Mammalia',     'Primates',      'Hominidae',    'Homo'),
        'Bos taurus':           ('Mammalia',    'Chordata',   'Mammalia',     'Artiodactyla',  'Bovidae',      'Bos'),
        'Gallus gallus':        ('Aves',        'Chordata',   'Aves',         'Galliformes',   'Phasianidae',  'Gallus'),
        'Amathia verticillata': ('Bryozoa',     'Bryozoa',    'Gymnolaemata', 'Ctenostomatida','Vesiculariidae','Amathia'),
        'Owenia fusiformis':    ('Annelida',    'Annelida',   'Polychaeta',   'Sabellida',     'Oweniidae',    'Owenia'),
    }
    
    def _make_df(species_list, site='Demo_Site'):
        rows = []
        for sp in species_list:
            dg, ph, cl, od, fm, gn = _base_tax[sp]
            rows.append({'Site': site, 'label': sp, 'display_group': dg,
                         'Phylum': ph, 'Class': cl, 'Order': od,
                         'Family': fm, 'Genus': gn,
                         'n_otus': len(species_list) + 5,
                         'n_otu_groups': 8})
        return pd.DataFrame(rows)
    
    # ---- Define which marker detects which species (edit freely) ----
    # These names become the legend labels and markers column values automatically
    marker_data = {
        'Vertebrates 12S': _make_df([
            'Dicentrarchus labrax', 'Sparus aurata', 'Diplodus vulgaris',
            'Sardina pilchardus', 'Chelon auratus', 'Mullus surmuletus',
            'Homo sapiens', 'Bos taurus', 'Gallus gallus',
            # Overlap with 18S:
            'Mytilus edulis', 'Carcinus maenas',
        ]),
        'Eukaryotes 18S': _make_df([
            'Patella caerulea', 'Patella vulgata', 'Octopus vulgaris',
            'Mytilus edulis', 'Carcinus maenas', 'Balanus trigonus',
            'Paracalanus parvus', 'Acartia clausii',
            'Paracentrotus lividus', 'Asterias rubens',
            'Amathia verticillata', 'Owenia fusiformis',
        ]),
        'Metazoa COI': _make_df([
            'Carcinus maenas', 'Balanus trigonus', 'Paracalanus parvus',
            'Mytilus edulis', 'Patella caerulea', 'Owenia fusiformis',
        ]),
    }
    
    # combine_marker_data builds 'markers' column dynamically from the dict keys
    combined_df, marker_names = combine_marker_data(marker_data)
    
    # Add metadata columns (normally comes from R/phyloseq)
    combined_df['Collector'] = 'J. Smith'
    combined_df['Location']  = 'Mediterranean Coast'
    combined_df['Date']      = '2024-06-15'
    
    fig, ax = create_wheel_of_life(
        combined_df[combined_df['Site'] == 'Demo_Site'],
        title='Demo Site (Multi-marker)',
        output_path='demo_wheel.pdf',
        marker_names=marker_names,
        collector='J. Smith',
        location='Mediterranean Coast',
        date='2024-06-15'
    )
    plt.close(fig)
    print("Done. Markers used:", marker_names)
    print("  (Change the marker_data dict keys above to use any names/count you like)")

