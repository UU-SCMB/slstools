from setuptools import setup, find_packages

setup(
    name="slstools",
    version="0.1.0",
    author="Roy Hoitink",
    author_email="L.D.Hoitink@uu.nl",
    license='GNU General Public License v3.0',
    long_description=open('README.md').read(),
    packages=find_packages(include=["slstools", "slstools.*"]),
    install_requires=[
        "numpy>=1.19.0",
        "miepython>=2.0.0",
        "scipy>=1.6.0",
    ],
)
