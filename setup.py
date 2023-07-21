import tinrun
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="tincheck",
    version=tinrun.VERSION,
    author="Aswathy Sebastian",
    author_email="aswathyseb@gmail.com",
    description="tincheck",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aswathyseb/tincheck",
    packages=find_packages(include=["tinrun", "tinrun.*"]),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],

    install_requires=[
        'plac',

    ],

    entry_points={
        'console_scripts': [
            'tincheck=tinrun.__main__:main',
        ],
    },

    python_requires='>=3.6',

)
