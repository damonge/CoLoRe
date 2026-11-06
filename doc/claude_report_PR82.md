# CoLoRe: faster line of sight for all the beam operations, and four bug fixes

Follow-up of PR 81 (`doc/claude_report_PR81.md`), on top of `master` (`7dc1ac8`):

1. **The overlapped ring is used for every beam operation**, not only for the source skewers (commit `3f4b2bb`).
2. **Only the fields that an active beam operation reads go around the ring** (commit `786c5ad`).
3. **Four bug fixes** found while testing (below): `write_cstm` (`bdfb6fe`), intensity maps without the Newtonian
   potential (`144cdb7`, a regression from PR 81), the velocity scaling of the skewers with `lpt_vels = 0` (`f6517a0`),
   and the lensing sums of the sources that were not set to zero (`d2890ce`).

The outputs of the beam operations other than the source skewers no longer match the previous ones bit by bit (the
sums get the slabs in a different order); they match to float rounding (below). No new parameters, no environment
variables.

## 1. Overlapped ring for all the beam operations

The line of sight sends the density (and velocity and Newtonian potential) slabs around a ring of the MPI tasks and
applies the beam operations to each slab. In PR 81 the ring in which the next slab is sent while the current one is
processed was used only for the source skewers, because for those every point gets contributions from at most two
slabs and the output stays bit-identical. Every other operation (kappa, lensing, ISW, intensity maps, custom maps)
adds the contribution of each slab to what it has accumulated, so the order of the slabs only changes the rounding of
the sums.

What changed (`src/beaming.c`):

- `ring_overlapped` takes the list of fields from the run: the density; the three velocities if `lpt_vels = 1`; the
  Newtonian potential if something reads it (`lpt_vels = 0`, kappa, ISW, lensing, or sources with lensing). The
  potential travels with its two neighbouring planes (`slice_left`, `slice_right`, which sit at the end of its
  buffer); the pointers to them are set again at every round.
- The per-slab calls of the operations are in one function (`beam_operations_one_slab`) used by both rings.
- The overlapped ring is used whenever there is more than one task. The original ring stays for a single task.

## 2. Bug in `write_cstm` (`src/io.c`)

Each task holds `num_pix = npix / NNodes` pixels of a custom map, listed in `listpix` (HEALPix nested indices).
`write_cstm` looped over all `npix` pixels of the full map and added `data[ip]` and `nadd[ip]` without `listpix`:
it read past the end of both arrays (garbage, sometimes NaN, in pixels the task does not hold) and put the local
pixels at the wrong positions. `write_kappa`, `write_isw` and `write_imap` already use `listpix`. Found with
AddressSanitizer (heap-buffer-overflow at `io.c` line 883). With one task `num_pix = npix` and the loop is harmless.
The fix uses `listpix` as the other writers do.

## What was checked

A 256^3 box (z from 0.05 to 1.0, 4 and 3 tasks on one node) with every output on at once: two source populations
(one with `include_lensing`), an intensity map, a custom map, kappa, ISW and, in the build with
`-D_USE_FAST_LENSING`, the lensing shells (nside 32 to 64). Master with the `write_cstm` fix against the new ring, for
the default build and the fast-lensing build, `lpt_vels = 1` and `0`: 8 runs, 17 to 24 files each.

- No NaN or garbage anywhere. The largest difference of any column, relative to the rms of the column of the master
  output, is 2.8e-6 (float32 resolution is 1.2e-7; the maps are sums over up to a thousand shells and slabs).
- The skewer files of the source population without lensing are byte-identical; the lensed one and all the maps
  agree to rounding.
- Smaller boxes (64^3 and 128^3 with 2, 3 and 4 tasks) agree to 2e-6 as well.
- The first comparison, before the `write_cstm` fix, showed garbage in the custom map in the new build and in the
  master build alike (the master one happens to read zeros or garbage depending on the memory layout); with the fix
  both agree.

**Timing: 2048^3 on 8 nodes** (8 tasks x 128 threads, one task per node, same parameters as above with nside 64: kappa,
ISW, custom and intensity maps, two source populations, one with lensing; `lpt_vels = 1`). Master (with the
`write_cstm` fix) and the new ring run one after the other in one job (60566323), one run each:

| | Master | New ring |
|---|---|---|
| Line of sight | 24.5 s | 11.2 s (2.2 times faster) |
| Whole run | 84.4 s | 74.5 s |

The master log shows 5.4 s of MPI time in the first slab alone; in the new ring task 0 waited 5.5 s for transfers in
total. The processing per slab is about 0.1 s with maps of nside 64, so most of the gain is the removed copy and wait,
and little of it is transfer hidden behind computation; with larger maps more of the transfer is hidden. Outputs: 27
FITS files compared; 8 byte-identical, 19 equal to rounding (largest difference 7.2e-6 of the rms of the column), no
NaN. Single measurements.

## Known limits

- The speed-up depends on how much of the line-of-sight time is the exchange: with small maps (nside 32 or 64 here)
  the processing per slab is short and the gain is mostly the removal of the copy and the wait; with large maps the
  overlap hides more of the transfer.
- The overlapped ring needs four to five more slabs of memory per task (one set of buffers per field).
- A ring with the potential and all the velocities sends five fields per round; fields that no active operation reads
  (the density when only kappa is on) are still sent, as in the original ring.

## 3. Only the fields that are read go around the ring

At 2048^3 on 8 nodes a round of the ring took about 1.15 s for about 21 GB per node (roughly 19 GB/s) while the
processing of a slab took 0.1 s, so the ring is limited by the bytes sent. `ring_fields_needed()` chooses the fields
from the operations that are on:

| Field | Sent when |
|---|---|
| Density | custom maps, or sources with stored skewers |
| LPT velocities | `lpt_vels = 1` and sources or fast lensing |
| Potential (with its two neighbouring planes) | kappa, ISW, fast lensing, sources with lensing, or sources with `lpt_vels = 0` |

The pointers of the fields that are not sent are NULL while the slabs are processed, so a read of one of them stops
the run instead of using the wrong slab. The log prints the fields that go around the ring. The intensity maps read
nothing in the line of sight (they use the potential before it).

### Timing: master (PR 81, with the four fixes) against this branch

2048^3 on 8 nodes (8 tasks x 128 threads), `lpt_vels = 1`, maps at nside 64, one run each:

| Operations on | Fields before -> after | Line of sight | Whole run |
|---|---|---|---|
| Kappa + ISW | 5 -> 1 | 23.3 s -> 2.6 s (8.9 times) | 78.2 s -> 57.4 s |
| Custom map only | 4 -> 1 | 14.1 s -> 2.5 s (5.6 times) | 64.9 s -> 53.7 s |
| Sources with lensing + kappa | 5 -> 5 | 20.3 s -> 11.4 s (1.8 times) | 78.0 s -> 71.8 s |
| Everything on | 5 -> 5 | 20.2 s -> 11.2 s (1.8 times) | 80.1 s -> 73.6 s |
| Sources with skewers only | 4 -> 4 | 8.1 s -> 8.2 s (same) | 60.6 s -> 60.8 s |

Master already had the overlapped ring for the source skewers, hence no change in the last row (output byte-identical).
Kappa + ISW and the custom map send fewer fields; the rest of the gain is the overlapped ring. A second run of the
custom map and of everything on gave 15.2 s -> 2.5 s and 21.4 s -> 11.8 s.

### Checks of the field selection

- 128^3, 4 and 3 tasks, `lpt_vels` 1 and 0, both builds (default and `-D_USE_FAST_LENSING`): the ring of the previous
  commit against this one, for kappa + ISW, custom maps, sources with and without skewers, sources with lensing and
  kappa, and everything on. All files are byte-identical except the intensity maps, which differ from run to run
  (about 2e-6) in master too.
- 2048^3 on 8 nodes: master with the fixes against this branch, files compared with `cmp` and, where they differ,
  column by column: worst difference 7.7e-6 of the rms of a column (everything on, including the intensity maps),
  4e-7 for kappa, ISW and the custom map.
- AddressSanitizer runs of the final code (128^3, 4 tasks): no errors for everything on with `lpt_vels = 0` and for
  the intensity maps alone.

## 4. Other bugs found while testing

- **Intensity maps without lensing, kappa or ISW (regression of PR 81).** `get_rvel` in `imap.c` always reads the
  Newtonian potential for the velocity of the cells, but `need_newtonian_potential()` did not count the intensity maps,
  so with only intensity maps and `lpt_vels = 1` the potential was not allocated and the run stopped with a
  segmentation fault. The fix adds `do_imap` to the function. DESI source runs have no intensity maps.
- **`MAX` instead of `MIN` in the velocity scaling of the skewers** (`srcs.c`, `lpt_vels = 0` only). The loop ran to
  `MAX(r*idr+0.5, nr-1)`, which can pass the end of the skewer of a source at the far end of the box (heap overflow,
  `malloc(): corrupted top size`). The loop that fills `v_skw` uses `MIN`. The output does not change, because the
  entries beyond are zero.
- **`kappa`, `dra` and `ddec` of the sources not set to zero** (default build, sources with `include_lensing`). The line
  of sight adds the lensing of each slab to them, but only `dz_rsd`, `e1` and `e2` were set to zero before, and the
  sources are allocated with `malloc`. Two runs of the same binary gave different values, from sane ones to ~1e3 for
  `KAPPA`. They are now set to zero with the others.
- **`write_cstm`**: section 2 above.

After these fixes two runs of the same binary agree byte for byte at 128^3 (everything on, `lpt_vels` 0 and 1), apart
from the intensity maps, whose files differ at the 1e-6 level from run to run (not investigated; it also happens in
master).
