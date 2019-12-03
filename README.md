# DynamicChangeBlindness


## Install psychopy for python 3.5 (ref: https://www.psychopy.org/installation.html)

```
conda install numpy scipy matplotlib pandas pyopengl pillow lxml openpyxl xlrd configobj pyyaml gevent greenlet msgpack-python psutil pytables requests[security] cffi seaborn wxpython cython pyzmq pyserial
conda install -c conda-forge pyglet pysoundfile python-bidi moviepy pyosf
pip install zmq json-tricks pyparallel sounddevice pygame pysoundcard psychopy_ext psychopy

pip install PyQt5
pip install pyglet==1.2.4
conda install jupyter
```

if pygame install fails with `SDL.h not found` (Mac OS), run

`pip install pygame==1.9.2`

if getting errors when saving movies

` pip install imageio==2.4.1`

------- Deprecated 

Uses python 2.7

## Install packages:

Install psychopy (ref: https://www.psychopy.org/installation.html)
With Anaconda:

```
conda create -n psypy python=2.7
conda activate psypy
conda install numpy scipy matplotlib pandas pyopengl pillow lxml openpyxl xlrd configobj pyyaml gevent greenlet msgpack-python psutil pytables requests[security] cffi seaborn wxpython cython pyzmq pyserial
conda install -c conda-forge pyglet pysoundfile python-bidi moviepy pyosf
pip install zmq json-tricks pyparallel sounddevice pygame pysoundcard psychopy_ext 
pip install psychopy 
```

## In case of install error

`ip install psychopy --user python `

Caution: install numpy *before* installing scipy.

TODO

Try these instead of official install https://psychologyit.uconn.edu/2017/09/20/instructions-for-installing-psychopy/
Try reinstalling for python 3.5

## References
Inspired by this paper:
Matlab library used in original paper: http://psychtoolbox.org/docs/