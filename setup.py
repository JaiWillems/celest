import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="celest",
    version="0.2.0",
    author="Jai Willems",
    author_email="jai52h@hotmail.com",
    description="Satellite dynamics and mission planning library.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JaiWillems/Celest",
    license="BSD-3-Clause",
    packages=setuptools.find_packages(include=["celest", "celest.*"]),
    install_requires=[
        "jplephem",
        "julian",
        "numpy",
        "pandas",
        "scipy"
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering"
    ],
    include_package_data=True,
    package_data={'': ['data/*.bsp']},
    python_requires='>=3.8',
)
