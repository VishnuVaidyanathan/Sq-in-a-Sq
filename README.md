# Sq in a Sq

Notes, data and tools from a study of the "squares in a square" problem:
pack n unit squares into the smallest possible square.

This repo holds work in progress. **Nothing here is a new record.**
The aim so far has been to understand why the best known packings look the way they do.

## Where the record data comes from

The best known packings are listed on David Ellsworth's page, which continues Erich Friedman's list:
<https://kingbird.myphotos.cc/packing/squares_in_squares.html>

All record side lengths in `data/` were read from that page (November 2026).
The three "record layout" files in `layouts/` are **not** the original files.
They were read off the published drawings and squeezed again with the tools here,
so they match the records only to about 0.001. Credit for those packings belongs to the people named on that page.

## What is in here

| Folder | What it holds |
|---|---|
| `data/` | Best known side for n = 1 to 324, written as m² + x and k² − r, with waste figures. Also an older floor-versus-best table for n = 1 to 100. |
| `layouts/` | Tile positions for each layout studied: one line per tile, `tile, x, y, angle`. |
| `pictures/` | A numbered picture of every main run. |
| `editor/` | A one-page editor (`index.html`). Open it in a browser to drag, turn and squeeze tiles. |
| `scripts/` | The Python scripts used. See the note at the end. |

In the n = 155 layouts marked "reference numbering", tile numbers are the same from file to file, and match the pictures.

## Findings so far

These come from reading the 147 tilted records up to n = 324. They are observations, not proofs.

1. **Each block starts tilted and ends straight.**
   Between two perfect squares, the tilted records start right after the lower square and run without a gap. Then the plain grid takes over up to the next square.

2. **The boundary creeps.**
   Write n = k² − r. Tilting beats the grid when r > k for k = 3 to 10, when r ≥ k for k = 11 to 15, and when r ≥ k − 1 for k = 16 to 18. The first case that breaks the "one missing row needs the full grid" guess is n = 110.

3. **About one row of waste.**
   For every tilted record, (side² − n) divided by m lies between 0.96 and 1.18, where m² is the perfect square below n. The average falls slowly, from 1.10 at m = 8 to 1.04 at m = 17.

4. **Ladders.**
   Seven decimal parts repeat from block to block (.5355, .5714, .6569, .7071, .7782, .8229, .8467) and cover 58 of the 147 records. Each ladder gains tiles at a steady rate. The .7071 ladder follows n = a² + a + 3 + ⌊(a−1)√2⌋ with side a + 1 + 1/√2.

5. **Soft ladders.**
   Most other records fall into groups whose decimal nearly repeats (near .60, .88, .93, .957, and .97 to .99).

6. **The channel rule (n = 131 and n = 155).**
   The small tilted block in the corner sits in a channel between the left wall and a flat wall of straight tiles. Its widest row sets the box:

   side = (number of straight columns) + 4·cos θ + sin θ

   With the block held at a fixed angle, this matched the squeeze to four decimals from −22° to −30.2°. It gives the .957 decimal at θ ≈ 30.3°.

7. **n = 181 works differently.**
   Its small block is near −26.5° and does not obey the channel rule. There the block pushes into the wall and the wall tiles turn a little.

8. **Bent blocks.**
   In the n = 155 layouts the large tilted block is not at one angle. Its columns fan over about 3°, which leaves thin slits adding up to about one tile of area. Forcing the block flat made the box 0.03 to 0.06 bigger.

## What was tried on n = 155 (record 12.95845)

| Layout | Side |
|---|---|
| Record layout, re-squeezed here | 12.9587 |
| 131-style join (wedge tile, main block starting with 2) | 12.9608 |
| Main block first row cut to 2 | 12.9621 |
| Small block held straight, best angle (−30.2°) | 12.9614 |
| Small block at −25°, tile 27 entered the 28/41 joint | 13.0010 |
| Tile 27 tip held in the 28/41 joint, pushed in 0.01 to 0.07 | 12.9689 to 13.0171 |

None of these beat the record. The squeeze used here only tightens a given layout. It cannot rearrange blocks, so these results say little about layouts that were not tried.

## Limits

- The squeeze is accurate to about 0.001.
- Side lengths were taken from the page as shown, mostly to 5 decimals.
- The "k² − 2 always needs side k" result is treated as proven in the older table. A 2026 paper disputes the original proof.
- `scripts/touch.py` (tiles 13 and 14 pushed up to touch tile 1) was written but its run was stopped before finishing. It is untested.

## About the scripts

They are working scripts, kept as they were run. Many read intermediate files by fixed names and will not run as they stand. The useful parts are:

- `seg.py`: reads tile positions from a screenshot of a packing drawing.
- `opt.py`: the squeeze (penalty method, L-BFGS, then a bisection on the side).
- `rb.py`, `rigid.py`: squeeze with one block held as straight rows.
- `pin.py`: squeeze with the tip of one tile held in a joint.
- `force2.py`: finds which tiles carry load in a jammed layout.
- `g.py`: the ladder and boundary analysis.

They need Python with numpy, scipy, autograd and opencv.
