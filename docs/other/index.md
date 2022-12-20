---
layout: single
classes: wide
title: Other topics
permalink: /other/
usemathjax: true

sidebar:
  nav: "other"
---

## Analytical Function Parser

OSIRIS includes an analytical function parser, so that the user can specify parameters for the simulation in the form of an analytical expression that is evaluated at run time. For details on the capabilities and syntax of the is function parser please check the [function parser](/other/function_parser) section.

## Grid diagnostics

Diagnostics of grid quantities in OSIRIS share a common interface that allows the user to do several different types of grid diagnostics besides the simple dumping of the grid quantity to disk. The user may choose to do data reduction operations (such as spatial averaging and lineouts) and also to perform time averaging operations to filter out higher (time) frequency components. For a full description of the syntax and available diagnostic types see the [grid diagnostics](/other/grid_diagnostics) section.