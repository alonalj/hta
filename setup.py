import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hta",
    version="0.1.0",
    author="Alona Levy-Jurgenson",
    author_email="levyalona@gmail.com",
    description="Statistically assess the level of both spatial and global heterogeneity within a spatial sample.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/alonalj/hta",
    project_urls={
        "Bug Tracker": "https://github.com/alonalj/hta/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "hta_stats"},
    packages=setuptools.find_packages(where="hta_stats"),
    python_requires=">=3.8",
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'matplotlib'
    ]
)