---
title: Finite element seminar
---

This is the page of **Math 692 Graduate Seminar: Finite elements** (CRN 35370) in Spring 2024, in the [Dept. of Mathematics and Statistics](http://www.uaf.edu/dms/) at the [University of Alaska Fairbanks](http://www.uaf.edu/).

* Course details:
  * Organizer/Instructor: [Ed Bueler](http://bueler.github.io/), [elbueler@alaska.edu](mailto:elbueler@alaska.edu).
  * Time and place: Thursdays 3:30-4:30pm, Chapman 206.  Please email me for the Zoom link.  In-person is preferred if you are on campus!
  * Credits: 1.0, but **non-credit attendance is welcomed**.

* [Firedrake installation advice for the seminar.](installation)

* The [Github repo for this website](https://github.com/bueler/fe-seminar) has example codes in the `py/` directory.

* [A schedule of topics is at the bottom.](#schedule)  This schedule is subject to change!

## guiding principles

1. show up
2. try stuff

## content

My plan is to take the lead on this seminar, and do most of the lecturing at the start.  Then to lead workshop style, or to encourage student presentations, as the semester goes on.  My goal is to teach the practical finite element method, through actual computations using [Firedrake](https://www.firedrakeproject.org/), in the first 6 or so weeks or so, and then branch out.  I plan to mix in enough theory to make the code make sense.  The goal is definitely _not_ to prove systematically that the finite element method converges for certain problems.

The basic scalar problems I have in mind are the [Poisson equation](https://en.wikipedia.org/wiki/Poisson%27s_equation), some [p-Laplacian](https://en.wikipedia.org/wiki/P-Laplacian) equations, and the [advection-diffusion equation](https://en.wikipedia.org/wiki/Convection%E2%80%93diffusion_equation).  Then we can also look at time-dependent problems, vector problems including the [Stokes equations for fluids](https://en.wikipedia.org/wiki/Stokes_flow), and even inequality-constrained problems like the [obstacle problem](https://en.wikipedia.org/wiki/Obstacle_problem).

## <a id="schedule"></a> topics (schedule)

| Date   | Topic | Slides | Code in repo |
|--------|-------|--------|--------------|
| 18 Jan | laptop day: install Firedrake, start Poisson | [slides](slides/18jan.pdf) | [code](https://github.com/bueler/fe-seminar/tree/main/py/)
| 25 Jan | Poisson equation examples, Paraview | [slide](py/25jan/poisson.pdf) | [code](https://github.com/bueler/fe-seminar/tree/main/py/25jan)
|  1 Feb | verification, assembly, errors | [slides](slides/1feb.pdf) | [code](https://github.com/bueler/fe-seminar/tree/main/py/1feb)
|  8 Feb | meshes & elements: Gmsh, P_k, Neumann | [slides](slides/8feb.pdf) | [code](https://github.com/bueler/fe-seminar/tree/main/py/8feb)
| 15 Feb | Stokes equations for glaciers | [slides](https://github.com/bueler/stokes-ice-tutorial/blob/main/slides.pdf) | [code (external)](https://github.com/bueler/stokes-ice-tutorial)
| 22 Feb | CANCELLED |
| 29 Feb | time-stepping | | [code](https://github.com/bueler/fe-seminar/tree/main/py/29feb)
|  7 Mar | fast solvers for Poisson? |
| 14 Mar | _Spring Break_ |
| 21 Mar | obstacle problem? |
| 28 Mar | advection-diffusion problems and DG? |
|  4 Apr | _student demo/presentation_ |
| 11 Apr | _student demo/presentation_ |
| 18 Apr | _student demo/presentation_ |
| 25 Apr | _student demo/presentation_ |
