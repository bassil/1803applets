from matplotlib.widgets import Slider, CheckButtons # import the Slider widget

import numpy as np
import matplotlib.pyplot as plt

class Matrix_Vector:
	def __init__(self, bounds, A, v):
		self.bounds = bounds
		self.A = A
		self.v = v
		self.Av_x, self.Av_y = self.A.dot(self.v)

		# Define the plot
		self.fig, self.ax = plt.subplots()
		# Adjust the plot to make room for our sliders and buttons
		plt.subplots_adjust(left=0.25, bottom=0.25)
		# We set the boundaries of our plot
		plt.axis([-self.bounds, self.bounds, -self.bounds, self.bounds])

		# generate eigenvalues and eigenvectors
		self.eigenvalues, self.eigenvectors = np.linalg.eig(self.A)
		# eigenvector 1
		self.x_v1, self.y_v1 = self.eigenvectors[:, 0]
		# eigenvector 2
		self.x_v2, self.y_v2 = self.eigenvectors[:, 1]

		# plot's of eigenvectors
		self.e1, = plt.plot([self.x_v1*-50, self.x_v1*50], 
			[self.y_v1*-50, self.y_v1*50], visible=False, 
			color='g', linewidth=0.4)
		self.e2, = plt.plot([self.x_v2*-50, self.x_v2*50], 
			[self.y_v2*-50, self.y_v2*50], visible=False,
			color='g', linewidth=0.4)

		# We draw the lines x = 0 and y = 0
		self.ax.axhline(y=0, color = 'k')
		plt.xticks(np.arange(-self.bounds, self.bounds+1, step=1))
		self.ax.axvline(x=0, color = 'k')
		plt.yticks(np.arange(-self.bounds, self.bounds+1, step=1))
		self.ax.set_axisbelow(True)
		self.ax.grid()

		# We draw the input vector in blue,
		self.input_vector = plt.quiver(0, 0, self.v[0], self.v[1], color='b', 
			angles='xy', scale_units='xy', scale=1, width = 0.005)
		# and the output vector in red
		self.output_vector = plt.quiver(0, 0, self.Av_x, self.Av_y, color='r', 
			angles='xy', scale_units='xy', scale=1, width = 0.005)
		

		# Define the sliders 
		self.axcolor = 'lightgoldenrodyellow'
		self.ax_x = plt.axes([0.25, 0.16, 0.65, 0.02], facecolor=self.axcolor)
		self.ax_y = plt.axes([0.25, 0.13, 0.65, 0.02], facecolor=self.axcolor)

		self.slider_x = Slider(self.ax_x, r'$x$', -6, 6, valinit=self.v[0])
		self.slider_y = Slider(self.ax_y, r'$y$', -6, 6, valinit=self.v[1])	

		self.ax_a11 = plt.axes([0.25, 0.1, 0.65, 0.02], facecolor=self.axcolor)
		self.ax_a12 = plt.axes([0.25, 0.07, 0.65, 0.02], facecolor=self.axcolor)
		self.ax_a21 = plt.axes([0.25, 0.04, 0.65, 0.02], facecolor=self.axcolor)
		self.ax_a22 = plt.axes([0.25, 0.01, 0.65, 0.02], facecolor=self.axcolor)

		self.slider_a11 = Slider(self.ax_a11, r'$a_{11}$', -6, 6, 
			valinit=self.A[0][0])
		self.slider_a12 = Slider(self.ax_a12, r'$a_{12}$', -6, 6, 
			valinit=self.A[0][1])
		self.slider_a21 = Slider(self.ax_a21, r'$a_{21}$', -6, 6, 
			valinit=self.A[1][0])
		self.slider_a22 = Slider(self.ax_a22, r'$a_{22}$', -6, 6, 
			valinit=self.A[1][1])

		# Define the checkboxes
		self.ax_check = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=self.axcolor)
		self.ax_check.set_title('Show eigen:')
		self.check = CheckButtons(self.ax_check, ('values', 
			'vectors', 'lines'), (False, False, False))

		# These are event handlers - The event is a slider value changing
		# if the slider value changes, run the update function.
		self.slider_x.on_changed(self.update)
		self.slider_y.on_changed(self.update)
		self.slider_a11.on_changed(self.update)
		self.slider_a12.on_changed(self.update)
		self.slider_a21.on_changed(self.update)
		self.slider_a22.on_changed(self.update)

		# This event is the clicking of the checkbox associated with 
		# eigenlines. 
		# TODO: showing the values of the eigenvalues and eigenvectors.
		# TODO: plotting decompositions
		self.check.on_clicked(self.check_box)


	# Update the values of A, v, and the eigenvalues
	def update(self, val):
		# updated value of v = (x, y)
		v_update = np.array([self.slider_x.val, self.slider_y.val])
		# updated value of A = ((a11, a12), (a21, a22))
		A_update = np.array([[self.slider_a11.val, self.slider_a12.val], 
			[self.slider_a21.val, self.slider_a22.val]])
		# updated value of Av
		Av_update = A_update.dot(v_update)
		x_update, y_update = v_update
		Av_x_update, Av_y_update = Av_update
		# updated eigenvalues and vectors
		eval_update, evec_update = np.linalg.eig(A_update)
		x_v1, y_v1 = evec_update[:, 0]
		x_v2, y_v2 = evec_update[:, 1]
		# Update values of inputvector, outputvector, and eigenlines
		self.input_vector.set_UVC(x_update, y_update)
		self.output_vector.set_UVC(Av_x_update, Av_y_update)
		self.e1.set_xdata([x_v1*-50, x_v1*50])
		self.e1.set_ydata([y_v1*-50, y_v1*50])
		self.e2.set_xdata([x_v2*-50, x_v2*50])
		self.e2.set_ydata([y_v2*-50, y_v2*50])
		self.fig.canvas.draw_idle()

	# If checked, display the eigenlines
	def check_box(self, label):
		if label == 'lines':
			self.e1.set_visible(not self.e1.get_visible())
			self.e2.set_visible(not self.e2.get_visible())
		plt.draw()

	# Finally, we show the plot
	def show(self):
		plt.show()


