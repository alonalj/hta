import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hta", # Replace with your own username
    version="1.0",
    author="Alona Levy-Jurgenson",
    author_email="levyalona@gmail.com",
    description="Spatial heterogeneity using HTA",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="TODO_https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)