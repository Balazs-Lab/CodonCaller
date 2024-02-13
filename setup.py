import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CodonCaller",
    version="1.0.0",
    author="Adam Nitido",
    author_email="adamnitido@gmail.com",
    description="Codon Based Variant Caller for HIV Deep Sequencing Data",
    long_description=long_description,
    url="https://github.com/adamnitido/CodonCaller",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
