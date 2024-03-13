# Libraries
import copy
import numpy as np
from abc import ABC, abstractmethod

# Abstract class for dynamical systems
class DynamicalSystem(ABC):

    def __init__(self, x0, t0, f, jac, dt, **kwargs):
        '''
        Instantiation of a dynamical system.
            Parameters:
                x0 (numpy.ndarray): Initial condition.
                t0 (float): Initial time.
                f (function): function f of ẋ = f(x,t) or x_(n+1) = f(x_n).
                jac (function): jacobian of f with respect to x.
                dt (float): time interval between two time steps.
                kwargs (dict): Dictionary of parameters for f and jac.
        '''
        self.x0 = x0
        self.t0 = t0
        self.x = x0
        self.t = t0
        self.dim = len(x0)
        self.f = f
        self.jac = jac
        self.dt = dt
        self.kwargs = kwargs

    def copy(self):
        '''
        Copy a dynamical system.
            Returns:
                A copy of the dynamical system.
         '''
        return copy.deepcopy(self)

    
    @abstractmethod
    def next(self):
        '''
        Compute the state of the system after one time step.
        '''
        pass
    
    @abstractmethod
    def next_LTM(self, W):
        '''
        Compute the state of a deviation vector after one time step.
            Parameters:
                W (numpy.ndarray): Array of deviations vectors.
            Returns:
                res (numpy.ndarray): Array of deviations vectors at next time step.
        '''
        pass
    
    def forward(self, n_steps, keep_traj):
        '''
        Forward the system for n_steps.
            Parameters:
                n_steps (int): Number of simulation steps to do.
                keep_traj (bool): Return or not the system trajectory.
            Returns:
                traj (numpy.ndarray): Trajectory of the system of dimension (n_steps + 1,self.dim) if keep_traj.
        '''
        if (keep_traj):
            traj = np.zeros((n_steps + 1,self.dim))
            traj[0,:] = self.x
            for i in range(1, n_steps + 1):
                self.next()
                traj[i,:] = self.x
            return traj
        else:
            for _ in range(n_steps):
                self.next()
    
# Continuous dynamical system
class ContinuousDS(DynamicalSystem):

    def __init__(self, x0, t0, f, jac, dt, **kwargs):
        '''
        Instantiation of a dynamical system.
            Parameters:
                x0 (numpy.ndarray): Initial condition.
                t0 (float): Initial time.
                f (function): function f of ẋ = f(x,t) or x_(n+1) = f(x_n).
                jac (function): jacobian of f with respect to x.
                dt (float): time interval between two time steps.
                kwargs (dict): Dictionary of parameters for f and jac.
        '''
        super().__init__(x0, t0, f, jac, dt, **kwargs)
    
    def next(self):
        '''
        Compute the state of the system after one time step with RK4 method.
        '''
        k1 = self.f(self.x, self.t, **self.kwargs)
        k2 = self.f(self.x + (self.dt / 2.) * k1, self.t + (self.dt / 2.), **self.kwargs)
        k3 = self.f(self.x + (self.dt / 2.) * k2, self.t + (self.dt / 2.), **self.kwargs)
        k4 = self.f(self.x + self.dt * k3, self.t + self.dt, **self.kwargs)
        self.x = self.x + (self.dt / 6.) * (k1 + 2*k2 + 2*k3 + k4)
        self.t += self.dt
    
    def next_LTM(self, W):
        '''
        Compute the state of a deviation vector after one time step with RK4 method.
            Parameters:
                W (numpy.ndarray): Array of deviations vectors.
            Returns:
                res (numpy.ndarray): Array of deviations vectors at next time step
        '''
        jacobian = self.jac(self.x, self.t, **self.kwargs)
        k1 = jacobian @ W
        k2 = jacobian @ (W + (self.dt / 2.) * k1)
        k3 = jacobian @ (W + (self.dt / 2.) * k2)
        k4 = jacobian @ (W + self.dt * k3)
        res = W + (self.dt / 6.) * (k1 + 2*k2 + 2*k3 + k4)
        return res

# Discrete dynamical system
class DiscreteDS(DynamicalSystem):

    def __init__(self, x0, t0, f, jac, dt = 1, **kwargs):
        '''
        Instantiation of a dynamical system.
            Parameters:
                x0 (numpy.ndarray): Initial condition.
                t0 (float): Initial time.
                f (function): function f of ẋ = f(x,t) or x_(n+1) = f(x_n).
                jac (function): jacobian of f with respect to x.
                dt (float): time interval between two time steps.
                kwargs (dict): Dictionary of parameters for f and jac.
        '''
        super().__init__(x0, t0, f, jac, dt, **kwargs)
    
    def next(self):
        '''
        Compute the state of the system after one time step using f.
        '''
        self.x = self.f(self.x, self.t, **self.kwargs)
        self.t += self.dt

    def next_LTM(self, W):
        '''
        Compute the state of a deviation vector after one time step with RK4 method.
            Parameters:
                W (numpy.ndarray): Array of deviations vectors.
            Returns:
                res (numpy.ndarray): Array of deviations vectors at next time step
        '''
        jacobian = self.jac(self.x, self.t, **self.kwargs)
        res = jacobian @ W
        if (self.dim == 1):
            return np.array([res])
        else:
            return res
