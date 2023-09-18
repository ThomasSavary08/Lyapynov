# Lyapynov

Lyapynov is a Python library to compute Lyapunov exponents, covariant Lyapunov vectors (CLV) and their adjoints for a dynamical system.
The results/algorithms used are taken from [P. Kuptsov's paper on covariant Lyapunov vectors](https://arxiv.org/abs/1105.5228).

## Installation

Use the package manager pip to install Lyapynov.

```bash
pip install Lyapynov
```

## Usage

First, one needs to define the system to study (which can be discrete or continuous) using the ContinuousDS or DiscreteDS methods. These methods take the following parameters as input:
* the initial conditions $x_{0}$ and $t_{0}$ $(x(t_{0}) = x_{0})$.
* the function $f$ describing the dynamical system:
```math
\left\{
    \begin{array}{ll}
        \dot{x} = f(x,t) \\
        x_{n+1} = f(x_{n},n)
    \end{array}
\right. 
```
* the jacobian of $f$ with respect to $x$ or $x_{n}$:
$$ \left\{ \begin{array}{ll} J(x,t) = \displaystyle \frac{\partial f}{\partial x}(x,t) \\ ~ \\ J(x_{n},n) = \displaystyle \frac{\partial f}{\partial x_{n}}(x_{n},n) \end{array} \right. $$

</br>

```python
import lyapynov

# Continous dynamical system
continuous_system = lyapynov.ContinuousDS(x0, t0, f, jac, dt)

# Discrete dynamical system
discrete_system = lyapynov.DiscreteDS(x0, t0, f, jac)
```

</br>
</br>

Once the dynamic system has been defined, the following functions can be used:

* Maximum Lyapunov exponents (MLE) to compute the maximum Lyapunov exponent $\lambda_{1}$.
```python
def mLCE(system : DynamicalSystem, n_forward : int, n_compute : int, keep : bool):
    '''
    Compute the maximal 1-LCE.
        Parameters:
            system (DynamicalSystem): Dynamical system for which we want to compute the mLCE.
            n_forward (int): Number of steps before starting the mLCE computation. 
            n_compute (int): Number of steps to compute the mLCE, can be adjusted using keep_evolution.
            keep (bool): If True return a numpy array of dimension (n_compute,) containing the evolution of mLCE.
        Returns:
            mLCE (float): Maximum 1-LCE.
            history (numpy.ndarray): Evolution of mLCE during the computation.
    '''
```

</br>

* Lyapunov characteristic exponents (LCE) to compute the $p$ first Lyapunov exponents $(\lambda_{1}, \cdots, \lambda_{p})$.
```python
def LCE(system : DynamicalSystem, p : int, n_forward : int, n_compute : int, keep : bool):
    '''
    Compute LCE.
        Parameters:
            system (DynamicalSystem): Dynamical system for which we want to compute the LCE.
            p (int): Number of LCE to compute.
            n_forward (int): Number of steps before starting the LCE computation. 
            n_compute (int): Number of steps to compute the LCE, can be adjusted using keep_evolution.
            keep (bool): If True return a numpy array of dimension (n_compute,p) containing the evolution of LCE.
        Returns:
            LCE (numpy.ndarray): Lyapunov Charateristic Exponents.
            history (numpy.ndarray): Evolution of LCE during the computation.
    '''
```

</br>


* Covariant Lyapunov vectors (CLV) to compute covariant Lyapunov vectors $\Gamma(t) = [\gamma_{1}(t), \cdots, \gamma_{m}(t)]$.
```python
def CLV(system : DynamicalSystem, p : int, n_forward : int, n_A : int, n_B : int, n_C : int, traj : bool, check = False):
    '''
    Compute CLV.
        Parameters:
            system (DynamicalSystem): Dynamical system for which we want to compute the mLCE.
            p (int): Number of CLV to compute.
            n_forward (int): Number of steps before starting the CLV computation. 
            n_A (int): Number of steps for the orthogonal matrice Q to converge to BLV.
            n_B (int): Number of time steps for which Phi and R matrices are stored and for which CLV are computed.
            n_C (int): Number of steps for which R matrices are stored in order to converge A to A-. 
            traj (bool): If True return a numpy array of dimension (n_B,system.dim) containing system's trajectory at the times CLV are computed.
        Returns:
            CLV (List): List of numpy.array containing CLV computed during n_B time steps.
            history (numpy.ndarray): Trajectory of the system during the computation of CLV.
    '''
```

</br>


* Adjoint covariant vectors (ADJ) to compute the adjoints of CLV $\Theta(t) = [\theta_{1}(t), \cdots, \theta_{m}(t)]$ such that $\Gamma(t)^{T} \Theta(t) = D(t)$ .
```python
def ADJ(CLV : list):
    '''
    Compute adjoints vectors of CLV.
        Parameters:
            CLV (list): List of np.ndarray containing CLV at each time step: [CLV(t1), ...,CLV(tn)].
        Returns:
            ADJ (List): List of numpy.array containing adjoints of CLV at each time step (each column corresponds to an adjoint).
    '''
```

</br>
</br>
An example of using the package is given in the notebook Example.ipynb.