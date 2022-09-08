# Particle Tracking

## Introduction

The particle tracking diagnostic saves detailed information on the time
evolution of simulation particles, namely position, momenta and energy.
To do so for every particle in the simulation would require an
impossible amount of storage space and furthermore most particles in a
given simulation will be of little interest. The solution is therefore
to track only "interesting" particles.

This will, however, pose an additional problem: the user/code must
choose which are the relevant particles to follow. In Osiris this is
addressed by running the simulation once, and choosing particles from
the diagnostic output. The user can then generate low resolution tracks
from the generated RAW particle diagnostics, or run the simulation again
providing Osiris with a list of particles to track.

## First Run

In order to track particles Osiris will have to attach a unique tag to
each particle in the simulation. This has a negligible impact on
performance (sligth communication/sorting overhead), but does require
extra memory (more 16 bytes per particle) so it is off by default. To be
able to generate tracks the user will therefore have to turn particle
tagging on. This is done in the corresponding species section:

```text
species
{
  ! (...)

  add_tag = .true.,     ! add tags to particles
}
```

For more information on enabling particle tagging check the [reference guide](../reference/Species.md).

Tagging information is included in the RAW diagnostics, so make sure you
also turn this diagnostic on. You can save a significant amount of
storage space using the RAW particle selection parameters:

```text
diag_species 
{
  ndump_fac_raw = 1,           ! do a RAW diagnostic at every ndump_fac iterations
  raw_gamma_limit = 10.0,  ! only save particles with a gamma >= 10
}
```

For more information on particle RAW diagnostics check the [reference guide](../reference/Species_Diagnostics.md).

## Post Processing

Particle selection is done using the particle data mining routines in
the visualization infrastructure. Once you have selected your
"interesting" particles, you can save the tags of those particles to
disk. Here's an example of such a tag file:

```text
! particle tag list
! original file : RAW-sp1-000019.hdf
! generated on Fri Jun 22 16:48:00 2007
! number of tags
          10
! particle tag list
        8,    467587
        7,  14108447
        7,  14067219
        7,  14018283
        7,   6702414
        7,   6151149
        8,   5774388
        8,   5724332
        8,  12931126
        7,    566385
```

Note that these are plain text files, so you can add or remove tags
manually. This allows you, for example, to select 2 distinct sets of
particles and then combine the 2 tag lists (just remember to set the
number of tags in the top of the file). Any line beginning with an
exclamation mark is ignored.

At this stage you can also use the visualization routines to create
tracks from the save RAW files. You are however limited to the frequency
of RAW dumps that you chose, so you will have a poor time resolution

## Second Run

Now that you know exactly which particles you want to track, you need to
run your simulation again, this time with tracking diagnostics on. To do
so you must provide Osiris with the tags list that you just generated:

```text
diag_species 
{
  ndump_fac_tracks = 10,         ! flush tracking information to storage at every 10 * ndump_fac iterations
  niter_tracks = 5,              ! store tracking data at every 5 iterations
  file_tags = "electrons.tags",  ! use the tags list from the file "electrons.tags" 
}
```

Please note that the numbers "ndump_fac_tracks" and "niter_tracks" are
independent. The "ndump_fac_tracks" should be similar to other
diagnostic "ndump_fac\*" parameters, which will result in faster
execution and smaller output files.

If you want you can change a few things in the input deck. The same
rules defined when restarting simulations apply. You can, for example,
remove all other diagnostics to speed things up. Just remember that
add_tag = .true. must still be present in the corresponding species
section.

For more information on particle tracking diagnostics parameters check
the [reference guide](../reference/Species_Diagnostics.md).

## Repacking Tracks Files (optional)

If you generate tracks files with thousands of tracks the visualization
routines can be quite slow in reading them. This is in part related to
the way that the HDF library organizes the data inside the tracks file,
since the file structure needs to change dynamically as the simulation
progresses, and can be significantly improved by repacking the file. For
each track the complete data is read into memory, and saved contiguously
inside the new file:

```text
 h5repack -v -i in-tracks.h5 -o out-tracks.h5 -l CONTI
```

For \~ 300 MB tracks file, with \~ 3000 tracks, the repacking process
took 6m24s and it lowered visualization I/O time from 32 sec. to 11 sec.
For more information see the
[h5repack](http://hdf.ncsa.uiuc.edu/HDF5/doc/Tools.html#Tools-Repack)
reference on the [HDF group](http://hdf.ncsa.uiuc.edu/) website. If you
don't have access to the HDF5 tools you can also do this from within IDL
using the *tracks_repack* command:

```text
IDL> tracks_repack 
```

This will take slightly longer than the above command (\~10 minutes for
the same test file as above), and the end files are not as efficient as
those produced by the HDF5 tool, but has the advantage of creating
smaller final files.
