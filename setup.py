from setuptools import setup, find_packages

setup(
    name='Lyapynov',
    version='0.1',
    author='Thomas Savary',
    author_email='savarythomas2102@gmail.com',
    description='A python package to compute Lyapunov exponents, covariant Lyapunov vectors (CLV) and adjoints of a dynamical system.',
    long_description=open('README.md').read(),  # Lisez le contenu du fichier README.md
    long_description_content_type='text/markdown',
    url='https://lien_vers_votre_projet',
    packages=find_packages(where='src'),  # Spécifiez le répertoire source
    package_dir={'': 'src'},  # Répertoire racine des packages
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    keywords='votre, mots-clés, ici',
    install_requires=[
        # Liste des dépendances requises
    ],
)
