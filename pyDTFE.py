import numpy as np
import os
path_current = os.path.abspath(os.path.dirname(__file__))

import sys
sys.path.append(path_current+"/build")
import src_dtfe



class DTFE_3D:
    def __init__(self, pos, 
                 L:float , 
                 Nmesh:int      = None, 
                 periodic:bool  = True,
                 paddingRate:float = None, 
                 ):
        '''
        Parameters
        ----------
        pos : Ndarray, shape (Nparts, 3). The position of particles.
        L : float. The size of the simulation box.
        Nmesh : int. The number of uniform mesh grid.
        periodic : bool, optional. The periodic boundary condition of simulation box.
        paddingRate : float, optional. The padding rate for the boundary of the simulation box.
        '''
        if pos.shape[1] != 3:
            raise ValueError("Input Position & Velocity Array should be shape liked (Nparts, 3)")
        if pos.dtype != np.float32:
            pos = pos.astype(np.float32)
        self.pos = pos
        self.L = L
        self.Nmesh = Nmesh
        self.periodic = periodic
        if periodic:
            if paddingRate is None:
                self.paddingRate = 0.06
            else:
                self.paddingRate = paddingRate
        else:
            self.paddingRate = 0.0

    def GridVel(self, vel, nmesh=None):
        '''
        Interpolate the velocity field of particles to the 3D uniform mesh grid.
        Asume the particle mass is 1, and one can rescale the output density to the physical value.
        Parameters
        ----------
        vel : Ndarray, shape (Nparts, 3). The velocity of particles.
        nmesh : int, optional. The number of uniform mesh grid.
        '''
        if vel.dtype != np.float32:
            vel = vel.astype(np.float32)
        if vel.shape != self.pos.shape:
            raise ValueError("The shape of input velocity should be same as the shape of input position.")
        if nmesh is None:
            nmesh = self.Nmesh
        if nmesh is None:
            raise ValueError("Nmesh is not set.")
        return src_dtfe.DTFE_3D_Grid( self.pos, vel, nmesh, self.L, self.paddingRate )
    
    def SampleVel(self, vel, samplings:np.array ):
        '''
        Interpolate the velocity field of particles to the given sampling position.
        Parameters
        ----------
        vel : Ndarray, shape (Nparts, 3). The velocity of particles.
        samplings : int. The position of the interpolaton sampling point, shape as (Nsampling, 3).
        '''
        if vel.dtype != np.float32:
            vel = vel.astype(np.float32)
        return src_dtfe.DTFE_3D_SampleVel( self.pos, vel, samplings, self.L, self.paddingRate )

    def SampleScalar(self, scalar, samplings):
        '''
        Interpolate the scalar field value of particles to the given sampling position.
        Parameters
        ----------
        scalar : Ndarray, shape (Nparts). The scalar field value of particles.
        samplings : int. The position of the interpolaton sampling point, shape as (Nsampling, 3).
        '''
        if vel.dtype != np.float32:
            vel = vel.astype(np.float32)
        return src_dtfe.DTFE_3D_SampleScalar( self.pos, scalar, samplings, self.L, self.paddingRate )
    


class DTFE_2D:
    def __init__(self, pos, 
                 L:float , 
                 Nmesh:int      = None, 
                 periodic:bool  = True,
                 paddingRate:float = None, 
                 ):
        '''
        Parameters
        ----------
        pos : Ndarray, shape (Nparts, 2). The position of particles.
        L : float. The size of the simulation box.
        Nmesh : int. The number of uniform mesh grid.
        periodic : bool, optional. The periodic boundary condition of simulation box.
        paddingRate : float, optional. The padding rate for the boundary of the simulation box.
        '''
        if pos.shape[1] != 2:
            raise ValueError("Input Position & Velocity Array should be shape liked (Nparts, 2)")
        if pos.dtype != np.float32:
            pos = pos.astype(np.float32)
        self.pos = pos
        self.L = L
        self.Nmesh = Nmesh
        self.periodic = periodic
        if periodic:
            if paddingRate is None:
                self.paddingRate = 0.06
            else:
                self.paddingRate = paddingRate
        else:
            self.paddingRate = 0.0

    def GridVel(self, vel, nmesh=None):
        '''
        Interpolate the velocity field of particles to the 2D uniform mesh grid.
        Parameters
        ----------
        vel : Ndarray, shape (Nparts, 2). The velocity of particles.
        nmesh : int, optional. The number of uniform mesh grid.
        '''
        if vel.dtype != np.float32:
            vel = vel.astype(np.float32)
        if nmesh is None:
            nmesh = self.Nmesh
        if nmesh is None:
            raise ValueError("Nmesh is not set.")
        return src_dtfe.DTFE_2D_Grid( self.pos, vel, nmesh, self.L, self.paddingRate )


    