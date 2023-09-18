from setuptools import setup, find_packages

setup(
    name='lyapynov',
    version='0.1',
    author='Thomas Savary',
    author_email='savarythomas2102@gmail.com',
    description='A python package to compute Lyapunov exponents, covariant Lyapunov vectors (CLV) and adjoints of a dynamical system.',
    url='https://github.com/ThomasSavary08/Lyapynov',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    keywords='Lyapunov, Lyapunov exponents, LCE, Covariant Lyapunov vectors, CLV, Dynamical systems, ODE',
    install_requires=['numpy'],
)