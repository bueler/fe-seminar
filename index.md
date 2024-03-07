---
title: Finite element seminar
---

Welcome to **Math 692 Graduate Seminar: Finite Elements** in Spring 2024!  We are in the [Dept. of Mathematics and Statistics](http://www.uaf.edu/dms/) at the [University of Alaska Fairbanks](http://www.uaf.edu/).

* Course details:
  * Organizer/Instructor: [Ed Bueler](http://bueler.github.io/), [elbueler@alaska.edu](mailto:elbueler@alaska.edu).
  * Time and place: Thursdays 3:30-4:30pm, Chapman 206.  Please email me for the Zoom link.  In-person is preferred if you are on campus!
  * Credits (CRN 35370): 1.0, but **non-credit attendance is welcomed**.

* [Firedrake installation advice.](installation)

* The [Github repo for this website](https://github.com/bueler/fe-seminar) has example codes in the `py/` directory.

* [A schedule of topics is at the bottom.](#schedule)  This schedule is subject to change!

## guiding principles

1. show up
2. try stuff

## content

My plan is to take the lead on this seminar, and do lectures at the start.  Then to lead workshop style, encouraging student presentations, as much as possible.

We will cover the practical finite element method, through actual computations using [Firedrake](https://www.firedrakeproject.org/).  I plan to mix in enough theory to make the code make sense.  The goal is definitely _not_ to prove systematically that the finite element method converges for certain problems.

An introductory [partial differential equation](https://en.wikipedia.org/wiki/Partial_differential_equation) problem is (of course) the [Poisson equation](https://en.wikipedia.org/wiki/Poisson%27s_equation).  Then we can look at
  * time-dependent problems like the [heat equation](https://en.wikipedia.org/wiki/Heat_equation),
  * [advection-diffusion equations](https://en.wikipedia.org/wiki/Convection%E2%80%93diffusion_equation),
  * vector problems including the [Stokes equations for fluids](https://en.wikipedia.org/wiki/Stokes_flow),
  * nonlinear problems like the [p-Laplacian](https://en.wikipedia.org/wiki/P-Laplacian) equation,
  * inequality-constrained problems like the [obstacle problem](https://en.wikipedia.org/wiki/Obstacle_problem),
and any other models of interest to the participants.

## <a id="schedule"></a> topics (schedule)

| Date   | Topic | Slides | Code in repo |
|--------|-------|--------|--------------|
| 18 Jan | laptop day: install Firedrake, start Poisson | [slides](slides/18jan.pdf) | [code](https://github.com/bueler/fe-seminar/tree/main/py/)
| 25 Jan | Poisson equation examples, Paraview | [slide](py/25jan/poisson.pdf) | [code](https://github.com/bueler/fe-seminar/tree/main/py/25jan)
|  1 Feb | verification, assembly, errors | [slides](slides/1feb.pdf) | [code](https://github.com/bueler/fe-seminar/tree/main/py/1feb)
|  8 Feb | meshes & elements: Gmsh, P_k, Neumann | [slides](slides/8feb.pdf) | [code](https://github.com/bueler/fe-seminar/tree/main/py/8feb)
| 15 Feb | Stokes equations for glaciers | [slides](https://github.com/bueler/stokes-ice-tutorial/blob/main/slides.pdf) | [code (external)](https://github.com/bueler/stokes-ice-tutorial)
| 22 Feb | CANCELLED |
| 29 Feb | time-stepping | [slides](slides/29feb.pdf) | [code](https://github.com/bueler/fe-seminar/tree/main/py/29feb)
|  7 Mar | fast solvers for Poisson | | [code](https://github.com/bueler/fe-seminar/tree/main/py/7mar)
| 14 Mar | _Spring Break_ |
| 21 Mar | obstacle problem? |
| 28 Mar | advection-diffusion problems and DG? |
|  4 Apr | _student demo/presentation_ |
| 11 Apr | _student demo/presentation_ |
| 18 Apr | Michael Christoffersen: elasticity |
| 25 Apr | _student demo/presentation_ |
