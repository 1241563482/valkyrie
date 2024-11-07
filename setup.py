from setuptools import setup, find_packages
#from hotpp import __version__


with open('README.md') as f:
    long_description = f.read()


setup(
    name="valkyrie",
    #version=__version__,
    author="Yijie Zhu",
    author_email="zhuyijie@smail.nju.edu.cn",
    #url="https://git.nju.edu.cn/bigd4/tpn",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy",     
        "ase"
    ],
    #license="MIT",
    #description="MiaoNet: Moment tensor InterAggregate Operation Net",
    #long_description=long_description,
    #long_description_content_type="text/markdown",
    entry_points={"console_scripts": ["val = valkyrie.entrypoints.main:main"]},
)
