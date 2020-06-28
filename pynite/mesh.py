import numpy as np
def regular_mesh_1D(start,end,step):
    N=int((end-start)/step+1)
    mesh=np.linspace(start,end,N)
    return mesh
