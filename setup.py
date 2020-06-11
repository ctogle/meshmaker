import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="meshmaker",
    version="0.0.1",
    author="Curtis Ogle",
    author_email="curtis.t.ogle@gmail.com",
    description="procedural mesh generation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ctogle/meshmaker",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: MIT License",
    ],
    python_requires='>=3.6',
)
