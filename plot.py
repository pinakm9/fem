import matplotlib.pyplot as plt
import numpy as np
import subprocess, os, re, time
import matplotlib.tri as tri
import matplotlib as mpl
import matplotlib.cm as cm

class Plotter(object):
	 
	def __init__(self, mesh, id_ = ''):
		self.mesh = mesh
		self.x = self.mesh.coordinates()[:,0] 
		self.id = str(id_)
		self.cmd = ['ffmpeg', '-y', '-framerate', '15', '-i', 'file_type', '-c:v', 'libx264', '-r', '30', '-pix_fmt', 'yuv420p', 'out_file']
	 
	def plot(self, func, dir_, itr, time, xlabel = 'x', ylabel = 'y', show = False):
		if dir_[-1] != '/':
			dir_ = dir_ + '/'
		self.dir = dir_
		self.func = np.vectorize(func) # func must be vectorized
		self.y = self.func(self.x)
		plt.clf()
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		self.quantity = ylabel
		plt.plot(self.x, self.y, label = 'time = {:.4f}'.format(time))
		plt.legend()
		if show is True:
			plt.show()
		plt.savefig(self.dir + self.quantity + '_' + self.id + '_' + '{:04.0f}'.format(itr))
		print('Plot #{} generated'.format(itr), end = '\r')
		

	def purge(self):
	    for file in os.listdir(self.dir):
	        if re.search('^'+self.quantity + '_' + self.id + '_\d*.png$', file):
	            os.remove(os.path.join(self.dir, file))

	def create_video(self, delete = True):
		file_type = self.quantity + '_' + self.id + '_' + '%04d.png'
		out_file = self.quantity + '_' + self.id + '.mp4'
		self.cmd[5] = file_type
		self.cmd[-1] = out_file
		self.cwd = os.getcwd() + '/' + self.dir
		p = subprocess.Popen(self.cmd, stderr = subprocess.PIPE, stdout = subprocess.PIPE, cwd = self.cwd)
		time.sleep(3)
		p.kill()
		if delete is True:
			self.purge()

class Plotter2(Plotter):

	def __init__(self, mesh, id_ = ''):
		super().__init__(mesh, id_)
		self.y = self.mesh.coordinates()[:,1]

	def plot(self, func, dir_, itr, time, xlabel = 'x', ylabel = 'y', quantity = '', show = False):
		if dir_[-1] != '/':
			dir_ = dir_ + '/'
		self.dir = dir_
		self.func = np.vectorize(func) # func must be vectorized
		plt.clf()
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		self.quantity = quantity 
		p = plt.tricontourf(self.x, self.y, self.mesh.cells(), self.func(self.x, self.y), 20)
		plt.colorbar(p)
		plt.legend([quantity, 'time = {:.4f}'.format(time)])
		if show is True:
			plt.show()
		plt.savefig(self.dir + self.quantity + '_' + self.id + '_' + '{:04.0f}'.format(itr))
		print('Plot #{} generated'.format(itr), end = '\r')

	def stream(self, func, dir_, itr, time, xlabel = 'x', ylabel = 'y', quantity = '', show = False):
		if dir_[-1] != '/':
			dir_ = dir_ + '/'
		self.dir = dir_
		self.func = np.vectorize(func) # func must be vectorized
		plt.clf()
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)