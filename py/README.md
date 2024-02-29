# py/

This directory contains a getting-started Python program

  * `start.py`

and sub-directories for examples after the first week.

## getting-started demo

You will need to have Firedrake installed on your machine or running on Google colab.  See the "Installation Advice" part of [the seminar website](https://bueler.github.io/fe-seminar/) and/or the instructions on [the Firedrake download page](https://www.firedrakeproject.org/download.html).

Once you have Firedrake on your own machine you should activate the virtual environment

        $ source firedrake/bin/activate

Then you should be able to run the `start.py` program:

        (firedrake) $ python3 start.py

If everything is working then this will generate a figure window (close it!) and write a file `result.pvd`.  The file contents can be viewed with [Paraview](https://www.paraview.org/):

        (firedrake) $ paraview results.pvd

To clean up the directory do

        (firedrake) $ make clean

## DATE/ subdirectories

See the slides for the given week/date on the table at [the seminar website](https://bueler.github.io/fe-seminar/).
