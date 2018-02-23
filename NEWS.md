## bsplinePsd 0.2.0

Added an argument called k1 in the gibbs_bspline function.  This allows the user to specify the starting value for parameter k.  If well-chosen, this can speed up convergence significantly.  The default is set to 20, which works well on all of the cases I have come across.  If missing (NA), then a random integer between 5 and kmax will be selected as the starting value for k.
