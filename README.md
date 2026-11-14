# Sq in a Sq

Short notes from a study of the "squares in a square" problem:
pack n unit squares into the smallest possible square.

Nothing here is a new record. These are notes on what was looked at and what was learned, kept as a starting point for later work.

Record data came from David Ellsworth's page, which continues Erich Friedman's list (read November 2026):
<https://kingbird.myphotos.cc/packing/squares_in_squares.html>

## What we looked at

- The best known side for every n from 1 to 324.
- Each n written two ways: as m² + x (how far above the perfect square below) and as k² − r (how far below the next one).
- The drawings for n = 131, 155 and 181, read into tile positions and squeezed again.
- Many hand-made changes to the n = 155 packing, each squeezed and compared with the record (12.95845).

## What we found

These are observations from the data, not proofs.

1. **Each block starts tilted and ends straight.** Between two perfect squares, tilted packings start right after the lower square and run without a gap. Then the plain grid takes over.
2. **The boundary creeps.** Tilting beats the grid when r > k for k = 3 to 10, when r ≥ k for k = 11 to 15, and when r ≥ k − 1 for k = 16 to 18. The first case where a grid missing one full row is beaten is n = 110.
3. **About one row of waste.** For every tilted record, (side² − n) / m lies between 0.96 and 1.18. The average falls slowly as n grows.
4. **Ladders.** Seven decimal parts repeat from block to block: .5355, .5714, .6569, .7071, .7782, .8229, .8467. Each gains tiles at a steady rate. The .7071 ladder follows n = a² + a + 3 + ⌊(a−1)√2⌋ with side a + 1 + 1/√2.
5. **Soft ladders.** Most other records fall into groups whose decimal nearly repeats (near .60, .88, .93, .957, and .97 to .99).
6. **The channel rule (n = 131 and 155).** The small tilted block in the corner sits between the left wall and a flat wall of straight tiles. Its widest row sets the box: side = straight columns + 4·cos θ + sin θ. With the block held at a fixed angle this matched the squeeze to four decimals from −22° to −30.2°. It gives the .957 decimal at θ ≈ 30.3°.
7. **n = 181 works differently.** Its small block is near −26.5° and does not obey the channel rule. The block pushes into the wall and the wall tiles turn a little.
8. **Bent blocks.** In n = 155 the large tilted block is not at one angle. Its columns fan over about 3°, leaving thin slits that add up to about one tile of area. Forcing the block flat made the box 0.03 to 0.06 bigger.

## What we tried on n = 155

| Change | Side |
|---|---|
| Record layout, squeezed again by us | 12.9587 |
| Join copied from n = 131 (wedge tile, large block starting with 2) | 12.9608 |
| Large block's first row cut from 3 to 2 | 12.9621 |
| Small block held as straight rows, best angle (−30.2°) | 12.9614 |
| Small block at −25°, one tile pushed into the wall joint | 13.0010 |
| That tile's tip held in the joint and pushed in by 0.01 to 0.07 | 12.9689 to 13.0171 |

None of these beat the record.

## What we learned about method

- Reading pictures and patterns explains which kind of layout wins and why. It does not find the last decimals.
- A gap of 0.002 cannot be judged by eye. Only a long machine search can settle it.
- Our squeeze only tightens a given layout (accurate to about 0.001). It cannot rearrange blocks, so it says little about layouts we did not build.
- The record for n = 155 came from a simulated annealing search followed by exact polishing. Beating it needs the same kind of search, seeded with both layout kinds above (channel kind and pushed-into-the-wall kind).

## Open questions

- Why is one row of area the price of tilting, and what lets a packing go below it?
- Is there a formula behind the soft ladders, with the angle shifting a little from one n to the next?
- What stops the channel-kind small block from going steeper than about 30°?
- No full public table exists for n = 325 to 1000.
