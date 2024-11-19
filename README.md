# Vista - VISualization Tool for Analysis

This is a CFD post-processing code for TUD INCA

Author : Wencan WU @ TU Delft @ July 6 2022

## python env setting

```bash
conda create --name pp python=3.9.18      # python3
conda install -c conda-forge mpi4py       # mpi, or via pip

pip3 install numpy==1.24.3                # usually all versions works
pip3 install matplotlib==3.7.1            # ！！！ newer version cannot plot separate isolines
pip3 install pyarrow==17.0.0               # pandas 3.0 will require it.
pip3 install pandas==2.2.0
pip3 install scipy==1.10.1                # interpolation

pip3 install pytecplot==1.6.1              # Snapshot Class

pip3 install numpy-stl                     # geometry module needs stl

pip3 install vtk==9.2.2                    # vtk
pip3 install pyvista==0.44.1               # pyvista

pip3 install opencv-python                 # no requirement
pip3 install moviepy==1.0.3                # when converting png to mp4

pip3 install xvfbwrapper                   # for rendering figure with vtk without display

pip3 install sendgrid==6.11.0              # for sending email with python

# fortran module need to be compiled manually.

cd /path/to/vista/lib/
f2py -m form -c *.f90
```