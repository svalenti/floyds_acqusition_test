# floyds_acqusition_test

This scripts help understanding if the object was inside the slit

- blue square is the position of the object after that the astronemtry is solved
- green point is the nominal position of the slit where LCO will try to put the object
- red point is the slit position computed by me (almost the same alqorithm used by LCO)


1) you need to have installed ffmpeg to make the video

2) copy all these files in the directory with all the guiding images

3) copy floacqastrodef.py in your /lib/python2.7/site-packages/ 

4) run plotposition.py
   -  it will compute the slit position
   -  solve the astrometry (using floacq.py
   - make the png files
   - make the video
 
