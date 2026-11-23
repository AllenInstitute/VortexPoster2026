r"""Export everything `index.ipynb` needs, so the page renders without the store.

Run this once (and again whenever a new materialization lands). It reads the
parquet store and the poster's own modules, and writes small CSV/JSON files
into `data/` plus web-sized copies of the poster figures into `figures/`.

Why an export step rather than reading the store from the page
--------------------------------------------------------------
The page is for posterity. The store is 10 GB of parquet on one laptop, the
poster modules import `figlib`, `population` and `standard_transform`, and
`cached_arbors` wants a 600 MB skeleton cache. None of that can be committed
and none of it will still be on anyone's machine in three years. What *can*
be committed is a few hundred kilobytes of derived numbers, and that is all
the page actually plots: 1,355 soma positions, a 13 x 13 type matrix, two
flow tables and a density curve.

So this script is the only thing that needs the store. `index.ipynb` reads
`data/` and nothing else, which also means the page can be rebuilt by GitHub
Actions, where the store does not exist.

Usage
-----
    ..\.venv\Scripts\python.exe build_assets.py           # everything
    ..\.venv\Scripts\python.exe build_assets.py --data    # just data/
    ..\.venv\Scripts\python.exe build_assets.py --figures # just figures/

The four panels taken from the poster PDF (blocks 01, 07, 08) are NOT rebuilt
here: they were never produced by this repository, so there is no source to
re-plot them from. They are cropped from the final PDF by a one-off script and
committed; see `figures/README.md`.
"""
from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
POSTER = HERE.parent / "poster"
DATA = HERE / "data"
FIGURES = HERE / "figures"
sys.path.insert(0, str(POSTER))

DATA.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)


def _modules():
    """Imported lazily so `--figures` does not need the store at all."""
    import posterbase as pb                                       # noqa: F401
    import adjacency as adj
    import arbors as ab
    import figlib as fl
    import reach as rk
    import sankey as sk
    import somata as so
    return adj, ab, fl, rk, sk, so


# --------------------------------------------------------------------------
# data/
# --------------------------------------------------------------------------


def export_palette(ab, adj) -> None:
    """The poster's own colours and order, so the page cannot drift from it.

    Written out rather than re-derived in the notebook: `MESO_INHERITED` is
    built by splitting cell-type colours at import time, and a page that
    rebuilt it from a different `arbors` would quietly recolour the figures
    the static panels sit beside.
    """
    payload = {
        "label_colour": {k: str(v) for k, v in adj.LABEL_PALETTE.items()},
        "meso_colour": {k: str(v) for k, v in ab.MESO_INHERITED.items()},
        "label_order": list(adj.LABEL_ORDER),
        "excitatory_order": list(adj.EXCITATORY_ORDER),
        "family_order": list(ab.FAMILY_ORDER),
        "display_name": {k: ab.display_label(k) for k in adj.LABEL_ORDER},
    }
    (DATA / "palette.json").write_text(json.dumps(payload, indent=2),
                                       encoding="utf-8")
    print(f"  palette.json          {len(payload['label_colour'])} colours")


def export_somata(ab, fl, so, adj) -> pd.DataFrame:
    """One row per tracked column cell: position, depth and both vocabularies.

    Excitatory and inhibitory in one table with a `broad_type` column, so the
    page can draw them together or apart from one file. Positions are the
    layer-aligned micrometres the hero uses, not raw EM nanometres.
    """
    frames = []
    for broad in ("excitatory", "inhibitory"):
        frame = so.positioned(fl.FINAL_VERSION, broad_type=broad,
                              level="meso_type", require_label=True)
        frame = frame.assign(broad_type=broad)
        frames.append(frame)
    cells = pd.concat(frames)

    # `label` for the page is the matrix's vocabulary: meso_type for the
    # excitatory cells, block 04's families for the inhibitory ones. Keeping
    # the raw meso_type too, since the page offers both groupings.
    cells["group"] = adj.label_of(cells.meso_type)
    out = pd.DataFrame({
        "nucleus_id": cells.index,
        "broad_type": cells.broad_type.to_numpy(),
        "meso_type": cells.meso_type.astype(str).to_numpy(),
        "group": cells.group.astype(str).to_numpy(),
        "um_x": cells.um_x.round(1).to_numpy(),
        "um_z": cells.um_z.round(1).to_numpy(),
        "depth_um": cells.um_y.round(1).to_numpy(),
    })
    out = out[out.group != "nan"]
    out.to_csv(DATA / "somata.csv", index=False)
    print(f"  somata.csv            {len(out):,} cells "
          f"({out.broad_type.value_counts().to_dict()})")
    return out


def export_types(adj, fl) -> None:
    """The 13 x 13 type matrix, both normalisations, and the denominators."""
    cells = adj.roster(fl.FINAL_VERSION, verbose=False)
    edges = adj.edges(cells, fl.FINAL_VERSION, verbose=False)
    total = adj.type_matrix(cells, edges)
    out_share = adj.normalise(total, "row")
    in_share = adj.normalise(total, "column")
    share = adj.local_share(cells, total, fl.FINAL_VERSION)

    long = (total.stack().rename("synapses").reset_index()
            .rename(columns={"level_0": "pre", "level_1": "post"}))
    long["out_share"] = out_share.stack().to_numpy()
    long["in_share"] = in_share.stack().to_numpy()
    long["synapses"] = long.synapses.astype(int)
    long.to_csv(DATA / "type_matrix.csv", index=False)
    share.round(3).to_csv(DATA / "type_totals.csv", index_label="group")
    print(f"  type_matrix.csv       {len(long)} pairs, "
          f"{int(long.synapses.sum()):,} synapses")
    print(f"  type_totals.csv       {len(share)} groups")

    stats = adj.summary(cells, edges)
    stats = {k: (float(v) if isinstance(v, float) else int(v))
             for k, v in stats.items()}
    (DATA / "graph_stats.json").write_text(json.dumps(stats, indent=2),
                                           encoding="utf-8")
    print(f"  graph_stats.json      density {stats['density']:.4f}, "
          f"reciprocal x{stats['reciprocal_enrichment']:.1f}")


def export_flow(sk, fl) -> None:
    """Block 02's two proofreading-status flow tables, for a Sankey.

    Written long rather than as the wide matrix `flow_tables` returns. The
    index is the *baseline* status and the columns the final one, so an
    `index=False` write silently throws away which row is which — the first
    version of this file was four unlabelled rows of counts.
    """
    tables = sk.flow_tables(fl.BASELINE_VERSION, fl.FINAL_VERSION)
    for name, table in tables.items():
        safe = str(name).replace(" ", "_").replace("-", "_")
        long = (table.rename_axis(index="before", columns="after")
                .stack().rename("cells").reset_index())
        long = long[long.cells > 0]
        long.to_csv(DATA / f"flow_{safe}.csv", index=False)
        print(f"  flow_{safe}.csv".ljust(24)
              + f"{len(long)} flows, {int(long.cells.sum()):,} cells")


def export_reach(ab, adj, fl, rk) -> None:
    """Output-synapse density along the medial-lateral axis, per cell type.

    Read from the synapse table once and binned per type here, rather than
    calling `reach.synapse_positions` once per type: that function caches on
    the cell set, so eight types would re-read a 10.9-million-row parquet
    eight times for the same answer.
    """
    cells = adj.roster(fl.FINAL_VERSION, verbose=False)
    excitatory = cells[cells.broad_type == "excitatory"]
    root_type = pd.Series(excitatory.group.to_numpy()
                          if "group" in excitatory else
                          excitatory.label.to_numpy(),
                          index=excitatory.pt_root_id.to_numpy())

    table = pd.read_parquet(
        fl.STORE / f"v{fl.FINAL_VERSION}" / "synapses.parquet",
        columns=["pre_pt_root_id", "ctr_pt_position_x",
                 "ctr_pt_position_y", "ctr_pt_position_z"])
    mine = table[table.pre_pt_root_id.isin(set(root_type.index))]
    voxels = mine[["ctr_pt_position_x", "ctr_pt_position_y",
                   "ctr_pt_position_z"]].astype("float64").to_numpy()
    um = np.asarray(ab._layer_transform().apply(
        voxels * np.asarray(fl.config.VOXEL_NM, dtype="float64")))

    width = 25.0
    lo = np.floor(um[:, 0].min() / width) * width
    hi = np.ceil(um[:, 0].max() / width) * width
    edges = np.arange(lo, hi + width, width)
    centres = (edges[:-1] + edges[1:]) / 2.0
    labels = root_type.reindex(mine.pre_pt_root_id).to_numpy()

    rows = []
    for group in pd.unique(labels):
        counts, _ = np.histogram(um[labels == group, 0], bins=edges)
        rows.append(pd.DataFrame({"um_x": centres, "group": str(group),
                                  "synapses": counts,
                                  "per_um": counts / width}))
    density = pd.concat(rows, ignore_index=True)
    density = density[density.synapses > 0]
    density.round(3).to_csv(DATA / "reach_density.csv", index=False)
    print(f"  reach_density.csv     {len(density):,} rows, "
          f"{int(density.synapses.sum()):,} synapses, "
          f"{density.group.nunique()} types")

    boundary = float(rk.boundary_um())
    (DATA / "reach_meta.json").write_text(json.dumps(
        {"boundary_um": boundary, "bin_width_um": width,
         "synapses": int(density.synapses.sum()),
         "cells": int(len(excitatory))}, indent=2), encoding="utf-8")
    print(f"  reach_meta.json       V1/HVA boundary at {boundary:,.0f} µm")


# --------------------------------------------------------------------------
# figures/
# --------------------------------------------------------------------------

#: Poster PNGs the page shows, and the width to resample each to. The
#: originals run to 19 MB -- block 05a alone is 4,779 px of rasterised arbor --
#: which is right for a 56-inch board and wrong for a web page.
COPIES = {
    "block02/block02a_proofreading_flow.png": 1600,
    "block02/block02b_neuron_296503.png": 1400,
    "block02/block02c_hours_296503.png": 1400,
    "block02/block02d_output_gains_column.png": 1200,
    "block04/block04b_hero_row.png": 2400,
    "block04/block04c_hero_only.png": 1800,
    "block04/block04f_rogues_gallery.png": 2400,
    "block04/block04d_orient_xy.png": 500,
    "block04/block04e_orient_zy.png": 500,
    "block05/block05a_column_topdown.png": 2200,
    "block05/block05b_halo_reciprocal.png": 2400,
    "block05/block05c_coreg_added.png": 1500,
    "block05/block05d_orient_xz.png": 500,
    "block05/block05e_adjacency.png": 1500,
    "block05/block05f_dots.png": 2400,
    "block06/block06a_reach_fan.png": 2000,
    "block06/block06b_layer5_fan.png": 2000,
}

#: The exemplar renders, shown as a gallery the poster had no room for.
SOLO_WIDTH = 700

#: Panels this repository never produced, supplied as source images rather
#: than plotted. They were cropped out of the poster PDF at first, which was
#: a mistake twice over: the boxes were set by eye and clipped two of the
#: panels, and rendering a placed bitmap at 200 dpi only interpolates it — the
#: PDF holds no more detail than the original. These are the originals.
#: Overridable, because it lives outside the repository and a hard-coded
#: absolute path would both break for anyone else and publish one person's
#: folder layout. Set POSTER_SUPPLEMENT to point at it; the default is a
#: sibling of this folder, which is where it would sit in a clone.
SUPPLEMENT = Path(os.environ.get("POSTER_SUPPLEMENT",
                                 HERE.parent / "sfn_2026_supplement"))

#: source file -> (target name, width). The page's prose carries the captions
#: and the citation that the poster printed around these, so the bare figures
#: are what is wanted here.
SUPPLEMENT_COPIES = {
    "microns_cover_compressed.png": ("panel01_volume.png", 756),
    "oligo_and_microglia.png": ("panel07_glia.png", 2200),
}

#: `mitochondria_quantification.png` holds both of block 08's figures stacked
#: with a white gutter between them. The page shows them as two, with text
#: in between, so it is split — at a row found by looking for the widest blank
#: band in the middle half of the image, not at a fraction typed from a
#: thumbnail.
MITO_SOURCE = "mitochondria_quantification.png"
MITO_PARTS = (("panel08_organelle.png", 1800), ("panel08_mito.png", 2200))


def export_figures() -> None:
    from PIL import Image

    source = POSTER / "output"
    written = 0
    for relative, width in COPIES.items():
        path = source / relative
        if not path.exists():
            print(f"  MISSING {relative}")
            continue
        written += _resample(Image, path, FIGURES / Path(relative).name, width)
    for path in sorted((source / "block04").glob("block04g_solo_*.png")):
        if "uncropped" in path.name:
            continue
        written += _resample(Image, path, FIGURES / path.name, SOLO_WIDTH)
    written += _export_supplement(Image)
    print(f"  {written} figures written to figures/")


def _export_supplement(Image) -> int:
    """The panels supplied as images rather than plotted here."""
    if not SUPPLEMENT.exists():
        print(f"  SUPPLEMENT not found at {SUPPLEMENT}; leaving "
              f"panel01/panel07/panel08 as they are")
        return 0
    written = 0
    for name, (target, width) in SUPPLEMENT_COPIES.items():
        path = SUPPLEMENT / name
        if not path.exists():
            print(f"  MISSING {name}")
            continue
        written += _resample(Image, path, FIGURES / target, width)

    path = SUPPLEMENT / MITO_SOURCE
    if path.exists():
        written += _split_and_write(Image, path)
    else:
        print(f"  MISSING {MITO_SOURCE}")
    return written


def _split_and_write(Image, path: Path) -> int:
    """Cut the stacked block-08 figure at its own white gutter."""
    import numpy as np

    with Image.open(path) as image:
        flat = Image.new("RGB", image.size, "white")
        flat.paste(image, mask=image.getchannel("A")
                   if image.mode == "RGBA" else None)
        blank = ~(np.asarray(flat) < 246).any(axis=2).any(axis=1)

        # Widest run of blank rows in the middle half of the image.
        lo, hi = int(flat.height * 0.25), int(flat.height * 0.75)
        best, start, longest = None, None, 0
        for y in range(lo, hi + 1):
            if y <= hi and blank[y]:
                start = y if start is None else start
                continue
            if start is not None and y - start > longest:
                longest, best = y - start, (start, y)
            start = None
        if best is None:
            raise SystemExit(f"no gutter found in {path.name}")
        split = (best[0] + best[1]) // 2
        print(f"  {path.name}: gutter rows {best[0]}–{best[1]}, "
              f"splitting at {split}")

        written = 0
        for (target, width), box in zip(
                MITO_PARTS,
                [(0, 0, flat.width, split),
                 (0, split, flat.width, flat.height)]):
            part = flat.crop(box)
            if part.width > width:
                height = round(part.height * width / part.width)
                part = part.resize((width, height), Image.LANCZOS)
            part.save(FIGURES / target, optimize=True)
            print(f"  {target:38s} {part.width:5d} px  "
                  f"{(FIGURES / target).stat().st_size / 1e6:5.2f} MB")
            written += 1
    return written


def _resample(Image, src: Path, dst: Path, width: int) -> int:
    with Image.open(src) as image:
        if image.width > width:
            height = round(image.height * width / image.width)
            image = image.convert("RGBA").resize((width, height),
                                                 Image.LANCZOS)
        # Flattened onto white: the poster panels are saved on white anyway,
        # and an RGBA PNG of a 2,400 px figure is about twice the size.
        flat = Image.new("RGB", image.size, "white")
        flat.paste(image, mask=image.getchannel("A")
                   if image.mode == "RGBA" else None)
        flat.save(dst, optimize=True)
    print(f"  {dst.name:38s} {flat.width:5d} px  "
          f"{dst.stat().st_size / 1e6:5.2f} MB")
    return 1


# --------------------------------------------------------------------------


def main(argv: list[str]) -> None:
    want_data = "--figures" not in argv
    want_figures = "--data" not in argv

    if want_data:
        print("data/")
        adj, ab, fl, rk, sk, so = _modules()
        export_palette(ab, adj)
        export_somata(ab, fl, so, adj)
        export_types(adj, fl)
        export_flow(sk, fl)
        export_reach(ab, adj, fl, rk)
    if want_figures:
        print("figures/")
        export_figures()
    print("\ndone. `quarto render` next.")


if __name__ == "__main__":
    main(sys.argv[1:])
