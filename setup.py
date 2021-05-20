import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Celest", # Replace with your own username
    version="0.1.6",
    author="Jai Willems",
    author_email="jai52h@hotmail.com",
    description="Satellite orbital position representation library.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JaiWillems/Celest",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    package_data={'': ['data/*.bsp']},
    python_requires='>=3.6',
)
