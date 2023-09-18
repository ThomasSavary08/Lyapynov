from setuptools import setup, find_packages

with open("README.md", "r", encoding = "utf-8") as file:
    desc = file.read()

setup(
    name='lyapynov',
    version='1.0',
    author='Thomas Savary',
    author_email='savarythomas2102@gmail.com',
    description='A python package to compute Lyapunov exponents, covariant Lyapunov vectors (CLV) and adjoints of a dynamical system.',
    long_description=desc,
    long_description_content_type="text/markdown",
    url='https://github.com/ThomasSavary08/Lyapynov',
    packages=find_packages(),
    keywords='Lyapunov, Lyapunov exponents, LCE, Covariant Lyapunov vectors, CLV, Dynamical systems, ODE',
    python_requires='>=3.8',
    install_requires=['numpy', 'matplotlib'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3", 
    ]
)