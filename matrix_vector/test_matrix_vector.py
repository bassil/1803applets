from matrix_vector import *

## Define your variables here:

# Matrix A:
a11 = 1
a12 = 1
a21 = 1
a22 = 0

A = np.array([[a11, a12],[a21, a22]])

# Vector v:
x = 1
y = 1

v = np.array([x,y])

## Define the plot interval [-bounds, bounds] 
#  (bounds>15 and the graphics start to look wonky)
bounds = 10

## Create a "matrix vector plot" object
p = Matrix_Vector(bounds, A, v)

## Show the plot
p.show()

##  Note, 
# 	the vector v is in blue, and
#   the vector Av is in red. 