from astropy.io import fits
import numpy as np

def read_fits_map(fits_image_filename):
	'''
	Reads FITS file, returns the header and the data filed
	'''
	with fits.open(fits_image_filename) as hdul:
		hdul.info()
		data = hdul[0].data
	
	return hdul, data


def get_convert():
	'''
	Returns a frequency converion matrix, nf x nc
	'''
# Freq          FF              Synch       

	conv = np.array([
	[ 3.566004, -1.755189],
	[ 3.226891, -1.463452],
	[ 2.930915, -1.220071],
	[ 2.671300, -1.016092],
	[ 2.442521, -0.844419],
	[ 2.240045, -0.699387],
	[ 2.060124, -0.576445],
	[ 1.899642, -0.471907],
	[ 1.755991, -0.382777],
	[ 1.626978, -0.306599],
	[ 1.510748, -0.241355],
	[ 1.405724, -0.185374],
	[ 1.310561, -0.137268],
	[ 1.224105, -0.095879],
	[ 1.145362, -0.060236],
	[ 1.073474, -0.029524],
	[ 1.007695, -0.003053],
	[ 0.947378, 0.019762],
	[ 0.891954, 0.039416],
	[ 0.840928, 0.056334],
	[ 0.793863, 0.070877],
	[ 0.750374, 0.083356],
	[ 0.710120, 0.094038],
	[ 0.672801, 0.103153],
	[ 0.638147, 0.110900],
	[ 0.605921, 0.117452],
	[ 0.575909, 0.122957],
	[ 0.547919, 0.127546],
	[ 0.521782, 0.131332],
	[ 0.497342, 0.134414],
	[ 0.474462, 0.136879],
	[ 0.453017, 0.138801]])
	return conv


def get_testconvert():
	'''
	Returns a simple frequency converion matrix, nf x nc
	'''
	conv = np.array([[1.,3],[2.,2.],[3.,1.]])

	return conv

def get_frequencies():
	'''
	Retunrs an array with the frequencies in MHz
	'''
	freqs = np.array([
	 90.533,
	 94.100,
	 97.667,
	101.234,
	104.801,
	108.367,
	111.934,
	115.501,
	119.068,
	122.635,
	126.201,
	129.768,
	133.335,
	136.902,
	140.469,
	144.036,
	147.602,
	151.169,
	154.736,
	158.303,
	161.870,
	165.436,
	169.003,
	172.570,
	176.137,
	179.704,
	183.270,
	186.837,
	190.404,
	193.971,
	197.538,
	201.104])
	return freqs

def get_testfrequencies():
	'''
	Retunrs an array with the frequencies in MHz
	'''
	freqs = np.array([100.,200.,300.])
	return freqs

#def get_fwhm():
# FWHM in arcmin
#	fwhm = np.full((32), 5.0) 
# Sigma in deg
#	sigma = fwhm/60.0/(2.0*np.sqrt(2.0*np.log(2.0)))
#	return fwhm, sigma

def get_fwhm(nf,s):
	'''
	Returns an array with FWHMs arcmin (all equal to the input s) and a corresponding Sigma in deg
	'''
# FWHM in arcmin
	fwhm = np.full((nf), s) 
# Sigma in deg
	sigma = fwhm/60.0/(2.0*np.sqrt(2.0*np.log(2.0)))
	return fwhm, sigma
		

def calcbins(nx, ny):
	'''
	Calculates the distance of each pixel from the centre of the patch,
	numerates distances in ascending order (this number is called bin number) 
	and returns the array with the bin numbers for each pixel. Several pixels can have the same
	bin number, they can be averaged by azav()
	'''
	icentre = nx/2
	jcentre = ny/2
	idsq = np.zeros(int((icentre+1)*(jcentre+1)))
	ipix = np.zeros(int((icentre+1)*(jcentre+1)), dtype=int)
	ibin = np.zeros(int((icentre+1)*(jcentre+1)), dtype=int)
	bins = np.zeros((nx,ny), dtype=int)

# create pixel number and distance-squared arrays
	n = -1
	for i in range(int(icentre+1)):
	    for j in range(int(jcentre+1)):
        	n = n+1 
        	ipix[n] = n
        	idsq[n] = (i - icentre)*(i-icentre) + (j - jcentre)*(j-jcentre)

	nlast = n

# sort idsq array into ascending order and impose same order on ipix
	arr1inds = idsq.argsort()
	idsq_sorted = idsq[arr1inds]
	ipix_sorted = ipix[arr1inds]

# create bin number array ibin
	bin = -1
	idsq_prev = -1000
	for i in range(nlast+1):
	    if idsq_prev != idsq_sorted[i]:
	        bin = bin+1
	    idsq_prev = idsq_sorted[i]
	    ibin[i] = bin

# populate 2D bins array with bin numbers from ibin w.r.t. the distance from the centre
# This is a modified version of F77 subroutine with indices shifted from [1:nx, 1:ny] to [0:nx-1, 0:ny-1]
	for n in range(nlast+1):
	    p=ipix_sorted[n]+1
	    i=int((p-1)/(jcentre+1) + 1)
	    j=int(p-(i-1)*(jcentre+1))
	    bin=ibin[n]
#	    print(n, p, i,j,bin)
	    bins[i-1,j-1]=bin
	    inew=int(2*(icentre+1)-i)
	    jnew=int(2*(jcentre+1)-j)
	    ilog= ((inew>=1) and (inew<=nx))
	    jlog= ((jnew>=1) and (jnew<=ny))
	    if (ilog): 
	        bins[inew-1,j-1]=bin
	    if (jlog):
	        bins[i-1,jnew-1]=bin
	    if (ilog and jlog):
	        bins[inew-1,jnew-1]=bin

	topbin=ibin[nlast]

	return bins, topbin


def makebeam(nx, ny, sigma, skycell):
	'''
	Creates 2D Gaussian profile
	'''
	centrex = nx/2
	centrey = ny/2
	nf = sigma.shape[0]
	beam = np.zeros((nf,nx,ny))
	beamarea = np.zeros((nf))

	for i in range(nx):
		for j in range(ny):
			dist=(skycell**2)*((i-centrex)**2+(j-centrey)**2)
			beam[:,i,j]=np.exp(-dist/(2.*sigma**2))
			beamarea=beamarea+beam[:,i,j]
	beam=beam/beamarea[:,None,None]
	beamarea=beamarea*(np.pi*skycell/180.)**2.

	return beam, beamarea

def azav(array, nx, ny, maxbin, bins):
	'''
	Averages azimutally w.r.t. the central pixel
	Replaces 2D array with the averaged values.
	Returns 1D arrays with the averaged values w.r.t. to the bin number (summ)
	as well as the number of times that particular bin (sq. distance) is used for the averaging (numm).
	'''
	num =  np.zeros((maxbin+1), dtype=int)
	summ = np.zeros((maxbin+1))

	for i in range(nx):
		for j in range(ny):
			ibin = bins[i,j]
			num[ibin] += 1
			summ[ibin] += array[i,j]

	summ /=num

	for i in range(nx):
		for j in range(ny):
			ibin = bins[i,j]
			array[i,j]=summ[ibin]

	return summ, num


def make_cosbell(nx, ny):
	'''
	Calculate cos(x,y) 2D bell-shaped profile 
	'''
	ic=nx/2
	jc=ny/2
	wss = 0.0
	cosbell = np.zeros((nx,ny))
	for i in range(nx):
		for j in range(ny):
			fc=np.cos((np.pi/nx)*(i-ic)) * np.cos((np.pi/ny)*(j-jc))
			wss=wss+fc**2.
			cosbell[i,j] = fc

	wss=wss/(nx*ny)
	return cosbell, wss


def make_icf(maps, nx, ny, nc, maxbin, bins, usecross):
	'''
	Calculates an Intrinsic Correlation Function L from the component covariance matrix C
	using Cholesky decomposition, C=np.matmul(L,L.T)
	'''
	icf = np.zeros((nc,nc,maxbin+1))
	choff  = np.zeros((nc,nc))
	
	cosbell, wss = make_cosbell(nx,ny)
	maps = maps*cosbell[None,:]

	mapsfft = np.fft.fft2(maps)
	mapsfft = np.fft.fftshift(mapsfft, axes=(1,2))
	
	for l in range(nc):
		for m in range(l+1):
			print(l,m)
			wkspce=np.abs(mapsfft[l]*np.conjugate(mapsfft[m]))/wss
			summ, num = azav(wkspce, nx, ny, maxbin, bins)
			for ibin in range(maxbin+1):
				icf[l,m,ibin]=summ[ibin] #/2.
				icf[m,l,ibin]=icf[l,m,ibin]

	if usecross :
		for ibin in range(maxbin+1):
			try:
				choff = np.linalg.cholesky(icf[:,:,ibin])
			except:
				print("WARNING: choldc failed at mode = ",ibin)
				for l in range(nc):
					choff[l,:] = 0.0
					choff[l,l] = np.sqrt(icf[l,l,ibin])

			icf[:,:,ibin] = choff

	return icf


def vistohid(icf, bins, visfft):
	'''
	Transforms Fourier data into hidden variable domain, H = np.matmul(np.linalg.inv(L), V)
	'''
	hidfft = np.zeros(visfft.shape, dtype=complex)
	for i in range(visfft.shape[1]):
		for j in range(visfft.shape[2]):
			ibin = bins[i,j]
			icfmat = icf[:,:,ibin]
			icfmat = np.linalg.inv(icfmat)
			hidfft[:,i,j] = np.matmul(icfmat, visfft[:,i,j])
	return hidfft


def hidtovis(icf, bins, hidfft):
	'''
	Transforms hidden Fourier data into visible domain, V = np.matmul(L, H)
	'''
	visfft = np.zeros(hidfft.shape, dtype=complex)
	for i in range(hidfft.shape[1]):
		for j in range(hidfft.shape[2]):
			ibin = bins[i,j]
			icfmat = icf[:,:,ibin]
			visfft[:,i,j] = np.matmul(icfmat, hidfft[:,i,j])
	return visfft
	

def updaterl(maxbin, convert, icf):
	'''
	Create a matrix product of FCM and L, RL
	'''
	nf = convert.shape[0]	
	nc = convert.shape[1]
	rl = np.zeros((nf,nc,maxbin+1))
	for ibin in range(maxbin+1):
		rl[:,:, ibin]=np.matmul(convert,icf[:,:, ibin])

	return rl

def predict(skyfft,beamfft,rl,datapfft,bins,i,j):
	'''
	Updates predicted Fourier data for datapfft(:,i,j) pixel using hidden Fourier data point (skyfft), 
	beam window function (beamfft) and RL function (FCM*L)
	'''
	ibin = bins[i,j]
	datapfft[:,i,j]=beamfft[:,i,j]*np.matmul(rl[:,:,ibin],skyfft[:,i,j])

def predict_mode(skyfft_mode,beamfft_mode,rl_mode):
	'''
	Updates predicted Fourier data for a single pixel/mode in the Fourier space using hidden Fourier data point (skyfft_mode), 
	beam window function (beamfft_mode) and RL function (FCM*L) for this mode
	'''
	
	return beamfft_mode*np.matmul(rl_mode,skyfft_mode)
	

