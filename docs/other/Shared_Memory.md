---
layout: single
classes: wide
title: Shared Memory Parallelism
permalink: /other/shared_memory
usemathjax: true

sidebar:
  nav: "other"
---

## Introduction

In addition to the standard spatial decomposition, OSIRIS can use a hybrid distributed memory (MPI) / shared memory (OpenMP) algorithm that relies on the fact that, usually, the computer nodes used have several computing cores sharing the node memory. Inside one of these nodes, particles are divided evenly across all cores, which ensures a perfect load balance between them, and may improve performance in some scenarios. However, this algorithm incurs in a (usually negligible) overhead, so in some cases, this may hinder performance.

## Using Shared Memory Parallelism

The use of the shared memory algorithm is controlled through the *n_threads* parameter in the [node_conf](../reference/node_conf) section. The user must specify the 2 levels of parallelism being used:

1. the standard node decomposition, using the node_number parameter. The code will launch 1 process per node specified.
2. the number of threads that each process will split into.

The total number of cores used will (should) be the total number of nodes times the number of threads per node specified.

## Choosing the number of threads

The *n_threads* parameter should be between 1, which corresponds to the default distributed memory algorithm, and the total number of cores in a compute node. For example, on a BlueGene/P the maximum value for *n_threads* should be 4, on a Cray XT5 (Dual hex-core AMD Opteron 2435@2.6GHz) the maximum value should be 12.

The number of threads being used does not have to be the maximum possible. For example, inside a Cray XT5 node you can use (processes × threads) 12 × 1, 6 × 2, 4 × 3, 3 × 4, 2 × 6 and 1 × 12. Generally speaking, the user should choose the maximum number of threads that fit in a single CPU, so for the Cray XT5 this would be 2 processes × 6 threads inside a compute node, and for the BlueGene/P this would be 1 × 4.

However, the performance will depend on the particular problem and hardware. The overhead of the shared memory algorithm grows with the number of threads, and the particular hardware configuration (NUMA nodes) may further limit the maximum number of threads. On the other hand, increasing the number of threads will lower load imbalance, so your mileage may vary.

## Example

Here's an example of a node_conf section for a 3D run, using 576 nodes in the $x_1$ direction and 8 nodes in the $x_2$ and $x_3$ directions, and launching 6 threads per node. The run will use 221184 cores in total.

```text
node_conf  {
    ! 36864 distributed memory nodes total
    node_number(1:3) =  576, 8, 8,
    ! 6 threads per node
    n_threads = 6,
}
```
