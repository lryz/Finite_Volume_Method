import numpy as np
def regular_mesh_1D(start,end,step):
	"""
    Parameters
    ----------
    start: float 
    		interval first value
    end: float 
    		interval last value 
    step: float 
    		step lenght bewteen 2 values
    Returns
    -------
        list of (end-start)/step + 1 values between start and end
    """
	N=int((end-start)/step+1)
	mesh=np.linspace(start,end,N)
	return mesh
