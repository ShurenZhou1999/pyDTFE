import numpy as np
import sys
sys.path.append("./build")
import src_dtfe


class DTFE_3D:
    def __init__(self, pos, 
                 L:float , 
                 Nmesh:int      = None, 
                 periodic:bool  = True,
                 paddingRate:float = None, 
                 ):
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
        if vel.dtype != np.float32:
            vel = vel.astype(np.float32)
        if nmesh is None:
            nmesh = self.Nmesh
        if nmesh is None:
            raise ValueError("Nmesh is not set.")
        return src_dtfe.DTFE_3D_Grid( self.pos, vel, nmesh, self.L, self.paddingRate )
    
    def SampleVel(self, vel, samplings):
        if vel.dtype != np.float32:
            vel = vel.astype(np.float32)
        return src_dtfe.DTFE_3D_SampleVel( self.pos, vel, samplings, self.L, self.paddingRate )

    def SampleScalar(self, scalar, samplings):
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
        if vel.dtype != np.float32:
            vel = vel.astype(np.float32)
        if nmesh is None:
            nmesh = self.Nmesh
        if nmesh is None:
            raise ValueError("Nmesh is not set.")
        return src_dtfe.DTFE_2D_Grid( self.pos, vel, nmesh, self.L, self.paddingRate )


    