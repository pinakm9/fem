import matplotlib.pyplot as plt
import numpy as np
import subprocess, os, re, time
import matplotlib.tri as tri

class Plotter(object):
	 
	def __init__(self, mesh, left, right, xlabel = 'x', ylabel = 'y', id_ = ''):
		self.mesh = mesh
		self.x = self.mesh.coordinates()[:,0] 
		self.left = left
		self.right = right
		self.xlabel = xlabel
		self.ylabel = ylabel
		self.id = str(id_)
		self.cmd = ['ffmpeg', '-y', '-framerate', '15', '-i', 'file_type', '-c:v', 'libx264', '-r', '30', '-pix_fmt', 'yuv420p', 'out_file']
	 
	def plot(self, func, dir_, itr, time, show = False):
		if dir_[-1] != '/':
			dir_ = dir_ + '/'
		self.dir = dir_
		self.func = np.vectorize(func) # func must be vectorized
		self.y = self.func(self.x)
		plt.clf()
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)
		plt.plot(self.x, self.y, label = 'time = {:.3f}'.format(time))
		plt.legend()
		if show is True:
			plt.show()
		plt.savefig(self.dir + self.ylabel + '_' + self.id + '_' + '{:04.0f}'.format(itr))
		print('Plot #{} generated'.format(itr), end = '\r')
		

	def purge(self):
	    for file in os.listdir(self.dir):
	        if re.search('^'+self.ylabel + '_' + self.id + '_\d*.png$', file):
	            os.remove(os.path.join(self.dir, file))

	def create_video(self, delete = True):
		file_type = self.ylabel + '_' + self.id + '_' + '%04d.png'
		out_file = self.ylabel + '_' + self.id + '.mp4'
		self.cmd[5] = file_type
		self.cmd[-1] = out_file
		self.cwd = os.getcwd() + '/' + self.dir
		p = subprocess.Popen(self.cmd, stderr = subprocess.PIPE, stdout = subprocess.PIPE, cwd = self.cwd)
		time.sleep(3)
		p.kill()
		if delete is True:
			self.purge()
	

class Plotter2(Plotter):

	def __init__(self, mesh, left, right, xlabel = 'x', ylabel = 'y', id_ = ''):
		super().__init__(mesh, left, right, xlabel, ylabel, id_)
		self.y = self.mesh.coordinates()[:,1]
		self.triangulation = tri.Triangulation(self.x, self.y)

	def plot(self, func, dir_, itr, time, show = False):
		if dir_[-1] != '/':
			dir_ = dir_ + '/'
		self.dir = dir_
		self.func = np.vectorize(func) # func must be vectorized
		plt.clf()
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel) 
		z = self.func(self.x, self.y)
		p = plt.tricontourf(self.triangulation, z, 20)
		cbar = plt.colorbar(p)
		if show is True:
			plt.show()
		plt.savefig(self.dir + self.ylabel + '_' + self.id + '_' + '{:04.0f}'.format(itr))
		print('Plot #{} generated'.format(itr), end = '\r')

	def stream(self, func):
		if dir_[-1] != '/':
			dir_ = dir_ + '/'
		self.dir = dir_
		self.func = np.vectorize(func)