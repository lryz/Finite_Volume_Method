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

def reguar_mesh_2D(start_x,end_x,step_x,start_y,end_y,step_y):
    x = regular_mesh_1D(start_x,end_x,step_x)
    y = regular_mesh_1D(start_y,end_y,step_y)
    return x,y