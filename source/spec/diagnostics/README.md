# Tracks diagnostic file format

The current version of the tracks files, described below, was designed to minimize the number of disk accesses at write time, resulting in much higher I/O performance, as well as much smaller files. However, this makes reading in the track data by the visualization routines more complicated and memory demanding.

Version 2 files can be converted into version 1 files using the `tracks_repack` command in **visXD**. Currently this works by reading the entire tracks file into memory and then saving the data using version 1 format, which means that you may get into trouble for very large ( > 4 GB ) track files.

## Track format Version 2

The version 2 format stores all track data in two 2D datasets named `data`and `itermap` that are organized as follows:

* `data[ nquants, np_total ]` - Where `nquants` is the number of quantities stored in the tracks file, and `np_total` the total number of track points, for all particles being tracked, that are stored in the file.
* `iter_map[ 3, n_segs  ]` - Where `n_segs` is the total number of track segments stored in the file. This dataset specifies which points in the `data` dataset correspond to each track segment. For segment `i`, `iter_map[1,i]` is the track id of the track the segment belongs to,  `iter_map[2,i]` is the number of points in the segment, and   `iter_map[3,i]` the starting iteration for the segment.

The sum of `iter_map[ 2, 1:i-1 ]` gives the (0 based) start position of segment `i` data in the `data` dataset.

 