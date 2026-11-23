Everything here is written by `build_assets.py --figures`, from two places:

- `../poster/output/` — the figures this repository plots. Listed in `COPIES`.
- `sfn_2026_supplement/` beside the poster PDF — the panels it does not plot,
  supplied as source images. Listed in `SUPPLEMENT_COPIES`:
  `panel01_volume` (the MICrONS volume render), `panel07_glia` (the
  oligodendrocyte and microglia panel) and `panel08_organelle` /
  `panel08_mito`, which come from one stacked image split at its own white
  gutter.

One file is neither, and is committed as-is:

- `block04h_astrocytes.webp` — the 126-astrocyte column render, supplied as a
  file. WebP at 926 × 2000 and 0.22 MB; as PNG it is several times that. The
  page shows it at 330 px and the lightbox opens the file itself, so there is
  no second copy to keep in step.

### These were cropped out of the poster PDF until 20 November 2026

That was wrong twice over and is worth recording, because nothing in the
pipeline could catch either half. The crop boxes were read off a 50 dpi
overview by eye, and two were too tight — block 07 lost 30.8 pt from its left
edge, taking the left of the "Oligodendrocyte myelinating axon" title and the
edge of the myelin inset with it, and block 08 lost about 10 pt on the left
and 8 pt top and bottom. A crop that is 3% too small still looks like a
figure. And rendering a *placed bitmap* out of a PDF at 200 dpi only
interpolates it: the PDF held no more detail than the original, so the
2,600 px crops were enlargements of the same pixels.

The source images are 8,870 px wide where the crop was 2,669. They also carry
no page furniture — no section rule bleeding in at the top, no teal footer
band at the bottom, which the crops picked up at two edges and had to be
nudged away from.

What the sources do not carry is the text set around them on the board: block
01's caption and the Princeton / Allen Institute / Baylor logos with the
*Nature* citation. The page's own prose carries the caption and the citation
with its DOI; the logos are not reproduced.
