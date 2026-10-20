# CoLoRe: faster 2LPT boxes with less memory (branch `gpu-port`)

The branch is called `gpu-port` for historical reasons (the same name is used in the other repositories of this
work); **nothing here uses a GPU**. All the changes are CPU code (C, OpenMP and MPI). They are on top of
`fe5d40e` of `damonge/CoLoRe` `master`, and the output is byte-identical to that code in every comparison made so
far. There are no new parameters and no changes to the parameter file; a few optional environment variables are
listed below, and all of them are off by default.

## Result at production size

4096^3 on 24 nodes (one task of 128 threads per node), the template parameters of the DR2 FS mocks (`lpt_vels=1`),
same parameter file and seed, the original and the new code run one after the other in one job on the same machine:

| | Original (`fe5d40e`) | This branch |
|---|---|---|
| Total time of the run | 1044.6 s | 289.8 s (3.6 times faster) |
| Node-hours (24 nodes) | 7.0 | 1.9 |
| Peak memory per node (Slurm) | 181 GB | 136 GB |
| Output | | all 24 `out_srcs_s1_*.fits` files identical with `cmp` |

The 1018 s and 175 GB per node in the log of the production run agree with the original code here. An earlier run of
the final code took 281.9 s, with these stage times: Gaussian field 8.8 s (production log 16 s), 2LPT 163.8 s (402 s),
point sources 6.2 s, line of sight 84.0 s (568 s), writes 19.0 s. The original code does not print per-stage times,
so the stage comparison is against the production log and not the same machine. One run each; the times are single
measurements.

At 1024^3 on one node (4 tasks x 32 threads) the whole run took 119.0 s and 85.9 GB before, 84.6 s and 83.3 GB after.

## What changed

| Where | Change | Measured result |
|---|---|---|
| `density.c`, `pos_2_dens_2lpt` | The cloud-in-cell deposit of the 2LPT particles and velocities ran on one thread per task (one thread of the 128 in the production layout). The grid rows (z plane, y) are now shared among the threads; each thread adds only to its own rows and visits the particles in index order, so every cell is summed in the same order as before. | 28.7 s to 1.8 s at 1024^3; 2LPT stage 42.8 to 18.4 s |
| `beaming.c`, line of sight | The Newtonian potential is not sent around the ring when nothing reads it; the received slab is copied back with all the threads; the next slab is sent and received with non-blocking calls while the current one is processed (two sets of buffers are swapped, and the own slab goes first). Used for the source skewers only. | 94 to 39 s at 2048^3 on 8 nodes; 84 s at 4096^3 on 24 nodes |
| `fourier.c`, `common.c` | The Newtonian potential grid was always allocated, filled, Fourier transformed and normalised. It is now skipped when nothing reads it (`lpt_vels=1` and no lensing, kappa or ISW; function `need_newtonian_potential()`). The random numbers are drawn as before. | 2LPT-stage peak at 1024^3 65 to 48-52 GB; Gaussian stage 2.67 to 1.83 s |
| `density.c`, `lpt_2` | Frees the three particle-velocity arrays after the deposit (the original allocates nine arrays and frees six), and writes the particle positions and velocities over the displacement arrays they are computed from (six arrays instead of nine, packed afterwards in blocks of rows). | whole-run peak 96.3 to 83.3 GB; 2LPT-stage peak 55.8 to 48.5 GB at 1024^3 |
| `density.c`, `share_particles_2lpt` | The particles to send are selected and packed in parallel, in the same order as before; the serial scan of an empty boundary buffer is skipped. | positions to deposit 30.0 to 12.6 s at 2048^3 on 4 nodes |
| `common.c`, `common.h`, `srcs.c` | Skewer sizes are computed in `size_t`: see the first bug below. | |

## Two bugs in the original, fixed here

- **32-bit `int` overflow in `catalog_alloc`.** The skewer arrays are allocated with `calloc(nsrc*nr, ...)` where
  `nsrc` and `nr = n_grid/2` are `int`. With 12.3M sources the product overflows 2^31 when sources per task times
  `n_grid/2` exceeds 2.1e9: 2048^3 on 4 tasks (3.1e9) stopped with "out of memory" right after the source
  redistribution, and so did 1024^3 on 2 tasks. 4096^3 on 24 tasks (1.05e9) and 2048^3 on 8 tasks are not affected,
  but 4096^3 on fewer than 12 tasks would be. The sizes (and one similar offset in `srcs.c`) are now `size_t`.
- **A memory leak.** With the 2LPT velocities, the three particle-velocity arrays are allocated and never freed, so
  they stay in memory through the line of sight and the output (13% of the peak at 1024^3).

## Optional environment variables (all off by default)

| Variable | Effect |
|---|---|
| `COLORE_TIMING=1` | Prints elapsed-time stamps of the stages, per-round processing and MPI times of the line of sight, and flushes the log at once (a killed job otherwise loses it). |
| `COLORE_SEND_NPOT=1`, `COLORE_NO_OVERLAP=1`, `COLORE_SERIAL_COPY=1` | Switch back to the original line-of-sight exchange (send the potential, no overlap, one-core copy), for comparisons. |
| `COLORE_SHARE_BLOCK`, `COLORE_UNPAD_BLOCK` | Block sizes of the particle sharing and of the packing of the rows (defaults 2^22 particles and a 512 MB buffer). |

## What was checked

- The `out_srcs_s1_*.fits` files compared with `cmp` against the output of the unmodified code after each change:
  1024^3 on 4 tasks (4 files, also with the original exchange switched back on); 512^3 on 2 and 3 tasks (uneven
  slabs); 2048^3 on 8 nodes, original exchange against the new one (8 files); a full 2048^3 run on 4 nodes with the
  final code (4 files); 4096^3 on 24 nodes (24 files). All identical.
- The slabs reach each task in a different order in the new ring. Every skewer point gets contributions from at
  most two slabs (the two z planes of the trilinear interpolation) and a sum of two terms does not depend on the
  order. This holds only for the source skewers; the new ring is not used for kappa, lensing, ISW or intensity
  maps, which sum over many slabs.
- Settings of the tests: the template parameters of the DR2 FS mocks with `lpt_vels = 1`, and the Makefile default
  `-D_BIAS_MODEL_2`. Two further 512^3 comparisons on 4 nodes (original `fe5d40e` against this branch, 4 files each)
  gave identical output: one built with `-D_BIAS_MODEL_3` (the linear bias model that the DR2 FS install note
  uses) and `lpt_vels = 1`, and one with `lpt_vels = 0` (velocities from the Newtonian potential).
- Not done: a comparison with the production box 0 (parameter file `param_box-0.cfg`, seed 202600000, 24 tasks x
  128 threads). That parameter file has no `lpt_vels` line, and both `fe5d40e` and the version of 2026-02-20
  (`ef2cf6b`, which asks for `lpt_vzty`) stop with "Couldn't read variable" on it, so production was built from a
  different version of the code than the ones tested here. A byte comparison needs the exact production commit and
  a build with `-D_BIAS_MODEL_3`.

## Known problems and limits

- 4096^3 on 12 nodes hung until the time limit, with repeated `cxil_map: write error` messages (the Slingshot
  interface failing to map memory), at about 204 GB per node, early in the run and before any line-of-sight output.
  My guess is that the large MPI buffers (23 GB per array at 12 tasks) exceed what the interface can map; I have not
  tried the original code at 12 nodes, so I do not know that it is new.
- The speed-up of the line of sight depends on `lpt_vels`. With `lpt_vels = 0` the Newtonian potential is needed
  for the velocities, so it is still sent around the ring and the overlapped ring (which needs `lpt_vels = 1`) is not
  used. At 512^3 on 4 nodes the output is identical and the time is the same (24.9 s and 25.0 s); the 4096^3 results
  above were measured with `lpt_vels = 1` only.
- The overlapped ring needs four more slabs of memory per node (46 GB at 4096^3 on 24 nodes).
- The FFTs are now the largest part of the 2LPT stage (about 60% of it). A GPU (cuFFTMp) version was considered and
  not started: it would need 13 to 25 GPU nodes for the memory, could save at most about 0.5 node-hours of the roughly 2 a box
  now takes, and cannot be byte-identical.

## Commits

`3d8fc6b` parallel 2LPT deposit; `319e35b` timing prints; `b412863` no potential in the ring; `6a9bd3d` parallel
copy of the received slab; `6e3d323` overlapped ring; `371cf6f` `COLORE_SEND_NPOT` and the other switches;
`e971aa1` `size_t` skewer sizes; `43d25d4` no potential when unused; `a51a10d` free the velocity arrays; `0f503a3`
positions and velocities in place; `14dabb9`, `f520047` parallel particle sharing; `6412f57` log flushed under
`COLORE_TIMING`.
