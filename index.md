---
title: Finite element seminar
---

This is the page of **Math 692 Graduate Seminar: Finite elements** (CRN 35370) in Spring 2024, in the [Dept. of Mathematics and Statistics](http://www.uaf.edu/dms/) at the [University of Alaska Fairbanks](http://www.uaf.edu/).

[Installation advice is below](#installation) and [a schedule of topics is at the bottom.](#schedule)

Organizer/Instructor: [Ed Bueler](http://bueler.github.io/), [elbueler@alaska.edu](mailto:elbueler@alaska.edu).

Time and place: Thursdays 3:30-4:30pm, Chapman 206.  Please email me for the Zoom link.  In-person is preferred if you are on campus!

Credits: 1.0, but **non-credit attendance is also welcomed**.

## guiding principles

1. show up
2. try stuff

## content

My plan is to take the lead on this seminar, and do most of the lecturing at the start.  Then to lead workshop style, or to encourage student presentations, as the semester goes on.  My goal is to teach the practical finite element method, through actual computations using [Firedrake](https://www.firedrakeproject.org/), in the first 6 or so weeks or so, and then branch out.  I plan to mix in enough theory to make the code make sense.  The goal is definitely _not_ to prove systematically that the finite element method converges for certain problems.

The basic scalar problems I have in mind are the [Poisson equation](https://en.wikipedia.org/wiki/Poisson%27s_equation), some [p-Laplacian](https://en.wikipedia.org/wiki/P-Laplacian) equations, and the [advection-diffusion equation](https://en.wikipedia.org/wiki/Convection%E2%80%93diffusion_equation).  Then we can also look at time-dependent problems, vector problems including the [Stokes equations for fluids](https://en.wikipedia.org/wiki/Stokes_flow), and even inequality-constrained problems like the [obstacle problem](https://en.wikipedia.org/wiki/Obstacle_problem).

## <a id="installation"></a> installation advice

<i>Thanks Stefano!</i>

### Firedrake on Google Colab

1. Navigate to [fem-on-colab.github.io](https://fem-on-colab.github.io/) and choose the [Packages tab](https://fem-on-colab.github.io/packages.html)
2. Scroll down to firedrake.  Copy the block of Python under "Real Mode".
3. Go to Google colab at [https://colab.research.google.com/](https://colab.research.google.com/)
4. In the pop-up window choose "+ New notebook".
4. Paste the block of Python into the first cell.
6. Hit shift-enter.  Stuff will be generated. After about a minute the spinny thing should stop.
7. Into the next cell type `from firedrake import *` and shift-enter.  This should give no warnings or errors.  Look for a green checkmark!  You have Firedrake!
8. As a further test and demo, copying my `start.py` script from [https://raw.githubusercontent.com/bueler/fe-seminar/main/py/start.py](https://raw.githubusercontent.com/bueler/fe-seminar/main/py/start.py) into a cell and shift-enter should cause no warnings or errors, and it should produce a nice figure!

### Firedrake on Linux/MacOS
 - For Linux and MacOS users simply follow the instructions outlined on [www.firedrakeproject.org/download.html](https://www.firedrakeproject.org/download.html)
 - It is highly recommended that, before you begin, you deactivate Anaconda. There are a few ways to go about this and they are outlined in the ['System anti-requirements' section of the firedrake download page](https://www.firedrakeproject.org/download.html#system-anti-requirements).
 - It is recommended that MacOS users use a Python distribution that is installed by Homebrew. To install Homebrew, visit [brew.sh](https://brew.sh/). Then you can install python using `brew install python3`.

### Firedrake on Windows
 - For now Windows support requires the use of Windows Subsystem for Linux (WSL)
 - https://github.com/firedrakeproject/firedrake/wiki/Installing-on-Windows-Subsystem-for-Linux seems to be the most up to date guide for getting this setup.

## <a id="schedule"></a> topics (speculative schedule)

| Date   | Topic | Slides |
|--------|-------|--------|
| 18 Jan | laptop day: install [Firedrake](https://www.firedrakeproject.org/) and solve Poisson | [slides](slides/18jan.pdf)
| 25 Jan | theory for the Poisson equation |
|  1 Feb | element types: P_k and quadrilaterals |
|  8 Feb | practical meshing |
| 15 Feb | performance for the Poisson equation |
| 22 Feb | advection-diffusion: CG is imperfect |
| 29 Feb | DG? |
|  7 Mar | Stokes equation |
| 14 Mar | _Spring Break_      |
| 21 Mar | obstacle problem |
| 28 Mar | _student demo/presentation_ |
|  4 Apr | _student demo/presentation_ |
| 11 Apr | _student demo/presentation_ |
| 18 Apr | _student demo/presentation_ |
| 25 Apr | _student demo/presentation_ |
