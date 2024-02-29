---
title: FE Seminar Installation Advice
---

# Installation of Firedrake

If the advice below is not enough, you might go to [Firedrake's installation debugging flowchart](https://www.firedrakeproject.org/download.html#debugging-install-problems).

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
 - [github.com/firedrakeproject/firedrake/wiki/Installing-on-Windows-Subsystem-for-Linux](https://github.com/firedrakeproject/firedrake/wiki/Installing-on-Windows-Subsystem-for-Linux) seems to be the most up to date guide for getting this setup.
