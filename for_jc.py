from setup import *
# use agg backend to prevent figure from showing
import matplotlib
matplotlib.use('qtagg')



# create a heavyside step function on the interval of x between -1 and 1
# and plot it
x = np.linspace(-2, 2, 100)
y = np.heaviside(x, 0.5)
# create a hyperbolic tangent function on the interval of x between -1 and 1
# and plot it
np.tanh(x)

plt.plot(x, y)
plt.show()
