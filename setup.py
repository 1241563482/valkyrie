from setuptools import setup, find_packages
from valkyrie import __version__


with open('README.md') as f:
    long_description = f.read()


setup(
    name="valkyrie",
    #version=__version__,
    author="Yijie Zhu",
    author_email="zhuyijie@smail.nju.edu.cn",
    url="https://github.com/1241563482/valkyrie",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy",     
        "ase",
        "matplotlib",
        "phonopy",
        "periodictable"
    ],
    license="MIT",
    description="Valkyrie: APEX LEGENDS NEVER DIE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={"console_scripts": ["val = valkyrie.entrypoints.main:main"]},
)

